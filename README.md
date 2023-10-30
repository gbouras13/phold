# phold
Phage Annotations with Protein Structures

* `phold run` Works only with CUDA GPU enabled for now.

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

# probably best to have torch cpu version with bioconda
# https://pytorch.org/get-started/previous-versions/
# pip install torch

# for my gridion
# need it choose the right version of torch compatible with your CUDA version
# Driver 11060 - v11.6
mamba install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia

# for HPC A100 
# Driver 11060 - v11.6
mamba install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia

# for macos (needed even though not used with remote)
mamba install pytorch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 -c pytorch

```

## Usage

```
Usage: phold run [OPTIONS]

  Runs phold

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -i, --input PATH         Path to input file in Genbank format  [required]
  -o, --output PATH        Output directory   [default: output_phold]
  -t, --threads INTEGER    Number of threads to use with Foldseek  [default:
                           1]
  -p, --prefix TEXT        Prefix for output files  [default: phold]
  -f, --force              Force overwrites the output directory
  --model_dir PATH         Path to save ProstT5_fp16 model to.
  -m, --model_name TEXT    Name of model: Rostlab/ProstT5_fp16.  [default:
                           Rostlab/ProstT5_fp16]
  -d, --database PATH      Path to foldseek PHROGs or ENVHOGs database.
                           [required]
  --database_name TEXT     Name of foldseek PHROGs or ENVHOGs database.
                           [default: all_phrogs]
  -e, --evalue TEXT        e value threshold for Foldseek  [default: 1e-3]
  -s, --sensitivity FLOAT  sensitivity parameter for foldseek  [default: 9.5]
  --batch_size INTEGER     batch size of ProstT5  [default: 1]
```


## Databases

* PHROG_ProstT5_Foldseek_db_updated - 440k non-redundant PHROG proteins
* ENVHOG_ProstT5_Foldseek_db - 2.2M representative ENVHOG proteins
    * need to make the ENVHOG mapping tsv work


## `run` if you have a CUDA GPU available
```
phold run -i tests/test_data/pharokka.gbk -o test_output -t 16 -f -d PHROG_ProstT5_Foldseek_db_updated/
```

## `remote` if you don't have a CUDA gpu (will take longer)
```
phold remote -i tests/test_data/pharokka.gbk -o test_output -t 16 -f -d PHROG_ProstT5_Foldseek_db_updated/
```


## TBD

* Add glue code for outputs (we can discuss what to include here)
* Make ENVHOGS DB work (with the mapping).
* Test different databases (PHROGs ENVHOGS). Tradeoff between size, speed and sensitivity. 
* Foldseek params (discuss with martin and milot).
* Try larger databases foldseek has already? Uniprot 50? Probably not worth it? Maybe swissprot?
* Get CPU working - will need help from Milot.
* Other ProstT5 models?