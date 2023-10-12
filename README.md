# phold
Phage Annotations with Protein Structures

* Works only with CUDA GPU for now.

## To create the env

```
git clone "https://github.com/gbouras13/phold.git"
cd phold
pip install -e .

# need foldseek
# will add the python deps into pyproject.toml later

mamba create -n prostt5 python=3.9 foldseek pip
conda activate prostt5
pip install transformers
pip install sentencepiece

# https://pytorch.org/get-started/previous-versions/
#pip install torch

# for my gridion
# need it choose the right version of torch compatible with your CUDA version

mamba install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia
```

## Tests

Required Susie's toy_foldseek_db for now as `-d`.

```
phold run -i tests/test_data/pharokka.gbk -o test_output -t 16 -f -d toy_foldseek_db/
```