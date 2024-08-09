#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Code adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/predict_3Di_encoderOnly.py

"""

# import dependencies

import copy
import re
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from datasets import Dataset
from loguru import logger
from torch.nn import CrossEntropyLoss
from torch.utils.data import DataLoader
from tqdm import tqdm
from transformers import (DataCollatorForTokenClassification, T5EncoderModel,
                          T5Tokenizer)
from transformers.modeling_outputs import TokenClassifierOutput
from transformers.models.t5.modeling_t5 import (T5Config, T5PreTrainedModel,
                                                T5Stack)
from transformers.utils.model_parallel_utils import (assert_device_map,
                                                     get_device_map)

from phold.features.predict_3Di import write_predictions
from phold.utils.constants import FINETUNE_DIR

# Modifies an existing transformer and introduce the LoRA layers


class LoRAConfig:
    def __init__(self):
        self.lora_rank = 4
        self.lora_init_scale = 0.01
        self.lora_modules = ".*SelfAttention|.*EncDecAttention"
        self.lora_layers = "q|k|v|o"
        self.trainable_param_names = ".*layer_norm.*|.*lora_[ab].*"
        self.lora_scaling_rank = 1
        # lora_modules and lora_layers are speicified with regular expressions
        # see https://www.w3schools.com/python/python_regex.asp for reference


class LoRALinear(nn.Module):
    def __init__(self, linear_layer, rank, scaling_rank, init_scale):
        super().__init__()
        self.in_features = linear_layer.in_features
        self.out_features = linear_layer.out_features
        self.rank = rank
        self.scaling_rank = scaling_rank
        self.weight = linear_layer.weight
        self.bias = linear_layer.bias
        if self.rank > 0:
            self.lora_a = nn.Parameter(
                torch.randn(rank, linear_layer.in_features) * init_scale
            )
            if init_scale < 0:
                self.lora_b = nn.Parameter(
                    torch.randn(linear_layer.out_features, rank) * init_scale
                )
            else:
                self.lora_b = nn.Parameter(torch.zeros(linear_layer.out_features, rank))
        if self.scaling_rank:
            self.multi_lora_a = nn.Parameter(
                torch.ones(self.scaling_rank, linear_layer.in_features)
                + torch.randn(self.scaling_rank, linear_layer.in_features) * init_scale
            )
            if init_scale < 0:
                self.multi_lora_b = nn.Parameter(
                    torch.ones(linear_layer.out_features, self.scaling_rank)
                    + torch.randn(linear_layer.out_features, self.scaling_rank)
                    * init_scale
                )
            else:
                self.multi_lora_b = nn.Parameter(
                    torch.ones(linear_layer.out_features, self.scaling_rank)
                )

    def forward(self, input):
        if self.scaling_rank == 1 and self.rank == 0:
            # parsimonious implementation for ia3 and lora scaling
            if self.multi_lora_a.requires_grad:
                hidden = F.linear(
                    (input * self.multi_lora_a.flatten()), self.weight, self.bias
                )
            else:
                hidden = F.linear(input, self.weight, self.bias)
            if self.multi_lora_b.requires_grad:
                hidden = hidden * self.multi_lora_b.flatten()
            return hidden
        else:
            # general implementation for lora (adding and scaling)
            weight = self.weight
            if self.scaling_rank:
                weight = (
                    weight
                    * torch.matmul(self.multi_lora_b, self.multi_lora_a)
                    / self.scaling_rank
                )
            if self.rank:
                weight = weight + torch.matmul(self.lora_b, self.lora_a) / self.rank

            # Convert weight to half (float16) precision to match unput
            # weight = weight.half()
            # input = input.half()

            # convert both to float
            input = input.to(torch.half)
            weight = weight.to(torch.half)

            return F.linear(input, weight, self.bias)

    def extra_repr(self):
        return (
            "in_features={}, out_features={}, bias={}, rank={}, scaling_rank={}".format(
                self.in_features,
                self.out_features,
                self.bias is not None,
                self.rank,
                self.scaling_rank,
            )
        )


def modify_with_lora(transformer, config):
    for m_name, module in dict(transformer.named_modules()).items():
        if re.fullmatch(config.lora_modules, m_name):
            for c_name, layer in dict(module.named_children()).items():
                if re.fullmatch(config.lora_layers, c_name):
                    assert isinstance(
                        layer, nn.Linear
                    ), f"LoRA can only be applied to torch.nn.Linear, but {layer} is {type(layer)}."
                    setattr(
                        module,
                        c_name,
                        LoRALinear(
                            layer,
                            config.lora_rank,
                            config.lora_scaling_rank,
                            config.lora_init_scale,
                        ),
                    )
    return transformer


class ClassConfig:
    def __init__(self, dropout=0.2, num_labels=3):
        self.dropout_rate = dropout
        self.num_labels = num_labels


class T5EncoderForTokenClassification(T5PreTrainedModel):
    def __init__(self, config: T5Config, class_config):
        super().__init__(config)
        self.num_labels = class_config.num_labels
        self.config = config

        self.shared = nn.Embedding(config.vocab_size, config.d_model)

        encoder_config = copy.deepcopy(config)
        encoder_config.use_cache = False
        encoder_config.is_encoder_decoder = False
        self.encoder = T5Stack(encoder_config, self.shared)

        self.dropout = nn.Dropout(class_config.dropout_rate)
        self.classifier = nn.Linear(config.hidden_size, class_config.num_labels)

        # Initialize weights and apply final processing
        self.post_init()

        # Model parallel
        self.model_parallel = False
        self.device_map = None

    def parallelize(self, device_map=None):
        self.device_map = (
            get_device_map(len(self.encoder.block), range(torch.cuda.device_count()))
            if device_map is None
            else device_map
        )
        assert_device_map(self.device_map, len(self.encoder.block))
        self.encoder.parallelize(self.device_map)
        self.classifier = self.classifier.to(self.encoder.first_device)
        self.model_parallel = True

    def deparallelize(self):
        self.encoder.deparallelize()
        self.encoder = self.encoder.to("cpu")
        self.model_parallel = False
        self.device_map = None
        torch.cuda.empty_cache()

    def get_input_embeddings(self):
        return self.shared

    def set_input_embeddings(self, new_embeddings):
        self.shared = new_embeddings
        self.encoder.set_input_embeddings(new_embeddings)

    def get_encoder(self):
        return self.encoder

    def _prune_heads(self, heads_to_prune):
        """
        Prunes heads of the model. heads_to_prune: dict of {layer_num: list of heads to prune in this layer} See base
        class PreTrainedModel
        """
        for layer, heads in heads_to_prune.items():
            self.encoder.layer[layer].attention.prune_heads(heads)

    def forward(
        self,
        input_ids=None,
        attention_mask=None,
        head_mask=None,
        inputs_embeds=None,
        labels=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
    ):
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        outputs = self.encoder(
            input_ids=input_ids,
            attention_mask=attention_mask,
            inputs_embeds=inputs_embeds,
            head_mask=head_mask,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        sequence_output = outputs[0]
        sequence_output = self.dropout(sequence_output)
        # Convert the tensor_half to torch.float32 before the operation
        sequence_output_fl = sequence_output.to(torch.float32)
        logits = self.classifier(sequence_output_fl)

        loss = None
        if labels is not None:
            loss_fct = CrossEntropyLoss()

            active_loss = attention_mask.view(-1) == 1
            active_logits = logits.view(-1, self.num_labels)

            active_labels = torch.where(
                active_loss, labels.view(-1), torch.tensor(-100).type_as(labels)
            )

            valid_logits = active_logits[active_labels != -100]
            valid_labels = active_labels[active_labels != -100]

            valid_labels = valid_labels.type(torch.LongTensor).to("cuda:0")

            loss = loss_fct(valid_logits, valid_labels)

        if not return_dict:
            output = (logits,) + outputs[2:]
            return ((loss,) + output) if loss is not None else output

        return TokenClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )


def PT5_classification_model(num_labels, model_dir, half_precision=False):
    # Load PT5 and tokenizer
    # possible to load the half preciion model (thanks to @pawel-rezo for pointing that out)

    # sets the device

    # Torch load will map back to device from state, which often is GPU:0.
    # to overcome, need to explicitly map to active device

    global device

    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        dev_name = "cuda:0"
    else:
        logger.error("Running phold with the finetuned ProstT5 model requires a GPU")

    # logger device only if the function is called
    logger.info("Using device: {}".format(dev_name))

    # load
    model_name = "Rostlab/ProstT5_fp16"
    logger.info(f"Loading T5 from: {model_dir}/{model_name}")
    logger.info(f"If {model_dir}/{model_name} is not found, it will be downloaded.")
    model = T5EncoderModel.from_pretrained(model_name, cache_dir=f"{model_dir}/").to(
        device
    )
    tokenizer = T5Tokenizer.from_pretrained(model_name, cache_dir=f"{model_dir}/")

    # Create new Classifier model with PT5 dimensions
    class_config = ClassConfig(num_labels=num_labels)
    class_model = T5EncoderForTokenClassification(model.config, class_config)

    # Set encoder and embedding weights to checkpoint weights
    class_model.shared = model.shared
    class_model.encoder = model.encoder

    # Delete the checkpoint model
    model = class_model
    del class_model

    # Print number of trainable parameters
    model_parameters = filter(lambda p: p.requires_grad, model.parameters())
    params = sum([np.prod(p.size()) for p in model_parameters])
    print("ProstT5_Classfier\nTrainable Parameter: " + str(params))

    # Add model modification lora
    config = LoRAConfig()

    # Add LoRA layers
    model = modify_with_lora(model, config)

    # Freeze Embeddings and Encoder (except LoRA)
    for param_name, param in model.shared.named_parameters():
        param.requires_grad = False
    for param_name, param in model.encoder.named_parameters():
        param.requires_grad = False

    for param_name, param in model.named_parameters():
        if re.fullmatch(config.trainable_param_names, param_name):
            param.requires_grad = True

    # Print trainable Parameter
    model_parameters = filter(lambda p: p.requires_grad, model.parameters())
    params = sum([np.prod(p.size()) for p in model_parameters])
    print("ProstT5_LoRA_Classfier\nTrainable Parameter: " + str(params) + "\n")

    return model, tokenizer


def load_model(filepath, model_dir, num_labels=1, mixed=False):
    # Creates a new PT5 model and loads the finetuned weights from a file

    # load a new model
    model, tokenizer = PT5_classification_model(
        num_labels=num_labels, model_dir=model_dir, half_precision=mixed
    )

    # Load the non-frozen parameters from the saved file
    non_frozen_params = torch.load(filepath)

    # Assign the non-frozen parameters to the corresponding parameters of the model
    for param_name, param in model.named_parameters():
        if param_name in non_frozen_params:
            param.data = non_frozen_params[param_name].data

    return tokenizer, model


# Dataset creation
def create_dataset(tokenizer, seqs):
    tokenized = tokenizer(seqs, max_length=1024, padding=False, truncation=True)
    dataset = Dataset.from_dict(tokenized)

    return dataset


def get_embeddings_finetune(
    cds_dict: Dict[str, Dict[str, Tuple[str, ...]]],
    model_dir: Path,
    output_3di: Path,
    max_batch: int = 100,
    finetuned_model_path: Optional[str] = None,
    proteins_flag: bool = False,
) -> bool:
    """
    Generate embeddings and predictions for protein sequences using finetuned PhrostT5 model.

    Args:
        cds_dict (Dict[str, Dict[str, Tuple[str, ...]]]): nested dictionary containing contig IDs, CDS IDs and corresponding protein sequences.
        model_dir (Path): Directory containing the pre-trained model.
        output_3di (Path): Path to the output 3Di file.
        max_batch (int, optional): Maximum batch size. Defaults to 100.
        finetuned_model_path (Optional[str], optional): Path to the finetuned model weights. Defaults to None.
        proteins_flag (bool, optional): Whether the sequences are proteins. Defaults to False.

    Returns:
        bool: True if embeddings and predictions are generated successfully.
    """

    predictions = {}
    if finetuned_model_path is None:
        finetuned_model_path = Path(FINETUNE_DIR) / "Phrostt5_finetuned.pth"
    else:
        finetuned_model_path = Path(finetuned_model_path)

    # Put both models to the same device
    tokenizer, model = load_model(
        finetuned_model_path, model_dir=model_dir, num_labels=20, mixed=False
    )

    global device

    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        dev_name = "cuda:0"
    else:
        logger.error(
            "Running phold with the finetuned PhrostT5 model requires an available GPU."
        )

    # logger device only if the function is called
    logger.info("Using device: {}".format(dev_name))

    model.to(device)

    logger.info("Beginning PhrostT5 Finetuned predictions")

    # loop over each record in the cds dict
    fail_ids = []

    for record_id, cds_records in cds_dict.items():
        # instantiate the nested dict
        predictions[record_id] = {}

        seq_record_dict = cds_dict[record_id]
        seq_dict = {}

        # gets the seq_dict with key for id and the translation
        for key, seq_feature in seq_record_dict.items():
            # get the protein seq for normal
            seq_dict[key] = seq_feature.qualifiers["translation"]

        # sort sequences by length to trigger OOM at the beginning
        seq_dict = dict(
            sorted(seq_dict.items(), key=lambda kv: len(kv[1][0]), reverse=True)
        )

        batch = list()

        for seq_idx, (pdb_id, seq) in enumerate(seq_dict.items(), 1):
            # replace non-standard AAs
            seq = seq.replace("U", "X").replace("Z", "X").replace("O", "X")
            seq_len = len(seq)
            seq = " ".join(list(seq))

            batch.append((pdb_id, seq, seq_len))

        # run them all

        pdb_ids, seqs, seq_lens = zip(*batch)

        # Create Dataset
        prediction_set = create_dataset(tokenizer, list(seqs))

        # Make compatible with torch DataLoader
        prediction_set = prediction_set.with_format("torch", device=device)

        # For token classification we need a data collator here to pad correctly
        data_collator = DataCollatorForTokenClassification(tokenizer)

        # Create a dataloader for the test dataset
        batch_size = max_batch
        prediction_dataloader = DataLoader(
            prediction_set,
            batch_size=batch_size,
            shuffle=False,
            collate_fn=data_collator,
        )

        # Put the model in evaluation mode
        model.eval()

        # Make predictions on the test dataset
        phrostt5_predictions = []

        with torch.no_grad():
            for batch in tqdm(prediction_dataloader):
                input_ids = batch["input_ids"].to(device)
                attention_mask = batch["attention_mask"].to(device)
                # Add batch results(logits) to predictions, we take the argmax here to get the predicted class
                phrostt5_predictions += (
                    model(input_ids, attention_mask=attention_mask)
                    .logits.argmax(dim=-1)
                    .tolist()
                )

        for batch_idx, identifier in enumerate(pdb_ids):
            s_len = seq_lens[batch_idx]

            # get the prediction and trim off the padded residues
            pred = phrostt5_predictions[batch_idx]
            pred = pred[0:s_len]

            predictions[record_id][identifier] = (pred, None, None)

    write_predictions(predictions, output_3di, proteins_flag)

    return True
