# phold
Phage Annotations with Protein Structures

Note: works by default with GPU if CUDA GPU enabled, but will also work CPU only (approx 1s/protein) with `--cpu`.

## To create the env

```
git clone "https://github.com/gbouras13/phold.git"
cd phold
pip install -e .

# need foldseek
# will add the python deps into pyproject.toml later

mamba create -n prostt5 python=3.9 foldseek pip
conda activate prostt5

# install the python deps

pip install transformers
pip install sentencepiece

# probably best to have torch cpu version with bioconda and make the gpu install an extra
# need it choose the right version of torch compatible with your CUDA version
# instructions at https://pytorch.org/get-started/previous-versions/
# pip install torch

# for gridion
# Driver 11060 - v11.6
mamba install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia

# for HPC A100 
# Driver 11060 - v11.6
mamba install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia

# for macos
mamba install pytorch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 -c pytorch

```

## Usage

```
Usage: phold [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  citation     Print the citation(s) for this tool
  compare      Runs phold compare (Foldseek)
  createdb     Creates phold compatible Foldseek db from AA FASTA and 3Di...
  createphrog  Creates phold compatible PHROG or ENVHOG foldseek db using...
  predict      Runs phold predict (ProstT5)
  proteins     Runs phold proteins (ProstT5 on a multiFASTA input)
  remote       Runs phold predict using foldseek API and compare locally
  run          Runs phold predict (ProstT5) and comapare (Foldseek)
```

* On large datasets, it is best to run `phold predict` (with GPU) followed by `phold compare` (CPU)

```
Usage: phold predict [OPTIONS]

  Runs phold predict (ProstT5)

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Path to input file in Genbank format  [required]
  -o, --output PATH      Output directory   [default: output_phold]
  -t, --threads INTEGER  Number of threads  [default: 1]
  -p, --prefix TEXT      Prefix for output files  [default: phold]
  -f, --force            Force overwrites the output directory
  --model_dir PATH       Path to save ProstT5_fp16 model to.
  -m, --model_name TEXT  Name of model: Rostlab/ProstT5_fp16.  [default:
                         Rostlab/ProstT5_fp16]
  --batch_size INTEGER   batch size for ProstT5. 1 is usually fastest.
                         [default: 1]
  --cpu                  Use cpus only.
  --omit_probs           Do not output 3Di probabilities from ProstT5
```

```
Usage: phold compare [OPTIONS]

  Runs phold compare (Foldseek)

Options:
  -h, --help                      Show this message and exit.
  -V, --version                   Show the version and exit.
  -i, --input PATH                Path to input file in Genbank format
                                  [required]
  --predictions_dir PATH          Path to directory with output3di.faa
  --pdb                           Use if you have pdbs for the input proteins
                                  (with AF2/Colabfold).
  --pdb_dir PATH                  Path to directory with pdbs
  --unrelaxed                     Use unrelaxed top rank pdb. By default,
                                  relaxed will be used (if found).
  -o, --output PATH               Output directory   [default: output_phold]
  -t, --threads INTEGER           Number of threads  [default: 1]
  -p, --prefix TEXT               Prefix for output files  [default: phold]
  -f, --force                     Force overwrites the output directory
  -d, --database PATH             Path to foldseek PHROGs or ENVHOGs database.
                                  [required]
  --database_name TEXT            Name of foldseek PHROGs or ENVHOGs database.
                                  [default: all_phrogs]
  -e, --evalue TEXT               e value threshold for Foldseek  [default:
                                  1e-3]
  -s, --sensitivity FLOAT         sensitivity parameter for foldseek
                                  [default: 9.5]
  --mode [tophit|topfunction|distribution]
                                  Mode to parse results.  [default: tophit]
```


## Databases

* PHROG_ProstT5_Foldseek_db_updated - 440k non-redundant PHROG proteins
* ENVHOG_ProstT5_Foldseek_db - 2.2M representative ENVHOG proteins
    * need to make the ENVHOG mapping tsv work


## `predict` if you have a CUDA GPU available
```
phold predict -i tests/test_data/pharokka.gbk -o test_output -t 1 -f 
```

```
phold comapre --predictions_dir test_output -o output_phold -t 1 -f -d PHROG_ProstT5_Foldseek_db/
```




## TBD

* LoRA on ProstT5.
* Split the predict/compare commands better.
* Better annotate PHROGs.
* Validate ProstT5 vs colabfold.
* Add the end glue code (to split out GFFs, GBKs etc).
* Allow for specifying structures more easily and remove `--unrelaxed`.
