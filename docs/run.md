# Running `phold`

`phold` is split into two overall modules regardless of your input - `predict` and `compare`. The reason `phold` has been split as such is to enable optimal resource usage.

`predict` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to translate protein amino acid sequences to the 3Di token alphabet used by [foldseek](https://github.com/steineggerlab/foldseek). This module is greatly accelerated if you have a GPU available and is recommended

`compare` runs [foldseek](https://github.com/steineggerlab/foldseek) to compare these ProstT5 predictions to the phold database (or alternatively, if you have provided pre-generated .pdb format protein structures for you proteins). This module does not require a GPU.

## Input 

Most subcommands of `phold` takes as their input an entry Genbank formatted file that contains the output of [`pharokka`](https://github.com/gbouras13/pharokka) for your phage or phage contigs. This will be called `pharokka.gbk` in your pharokka output.

Alternatively, `phold` will detect if the input is a FASTA contig/genome file as input. [`Pyrodigal-gv`](https://github.com/althonos/pyrodigal-gv]) will then be run to quickly predict the CDS and these will be annotated. However, neither tRNAs, tmRNA nor CRISPR repeats will be predicted (unlike in pharokka).

For `phold proteins-predict` and `proteins-compare`, the input will be a FASTA format file containing protein sequences in Amino Acids. These modes are useful for annotating bulk phage proteins. 

## Subcommands

### `predict`

To run phold predict



```
Usage: phold predict [OPTIONS]

  Uses ProstT5 to predict 3Di tokens - GPU recommended

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Path to input file in Genbank format  [required]
  -o, --output PATH      Output directory   [default: output_phold]
  -t, --threads INTEGER  Number of threads  [default: 1]
  -p, --prefix TEXT      Prefix for output files  [default: phold]
  -f, --force            Force overwrites the output directory
  -m, --model_dir PATH   Path to save ProstT5_fp16 model to.
  --model_name TEXT      Name of model: Rostlab/ProstT5_fp16.  [default:
                         Rostlab/ProstT5_fp16]
  --batch_size INTEGER   batch size for ProstT5. 1 is usually fastest.
                         [default: 1]
  --cpu                  Use cpus only.
  --omit_probs           Do not output 3Di probabilities from ProstT5
  --finetune             Finetune
  --finetune_path TEXT   Path to finetuned model weights
  ```


