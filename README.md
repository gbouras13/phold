# phold
Annotating phages with structural homology

`phold` is sensititve gene annotation tool for bacteriophage genomes and metagenomes using strucutal homology. 

`phold` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to translate protein amino acids to the 3Di token alphabet used by [foldseek](https://github.com/steineggerlab/foldseek). Foldseek is then used to search these against a database of 803k phage protein structures predicted using [Colabfold](https://github.com/sokrypton/ColabFold). Optionally, you can specify protein structures for your phage(s) instead of using ProstT5.

# Table of Contents

- [phold](#phold)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Output](#output)
  - [Usage](#usage)
  - [Databases](#databases)

# Installation

The only way to install `phold` is from source for now. 

Pypi and conda installations will be available soon. 

The only required non-Python dependency is `foldseek`. To install `phold` in a conda environment using mamba:

```
mamba create -n pholdENV pip foldseek 
conda activate pholdENV
git clone https://github.com/gbouras13/phold.git
cd phold 
pip install -e .
```

`phold` will run in a reasonable time for small datasets with CPU only (e.g. <5 minutes for a 50kbp phage).

 However, `phold predict` will run faster if a GPU is available, and is necessary for large metagenome datasets to run in a reasonable time. 

To utilise `phold` with GPU, a GPU compatible version of `pytorch` must be installed. 

If it is not automatically installed via pip/conda, please see [this link](https://pytorch.org) for more instructions on how to install `pytorch`. If you have an older version of CUDA installed, then you might find [this link useful](https://pytorch.org/get-started/previous-versions/).


# Quick Start

* `phold` takes a GenBank format file output from [pharokka](https://github.com/gbouras13/pharokka) as its input.
* It is most efficient to run `phold` in 2 steps for optimal resource usage.

1. Predict the 3Di sequences with ProstT5 using `phold predict`. This is massively accelerated if a GPU available.
* If you do not have a GPU available, add `--cpu`

```
phold predict -i tests/test_data/pharokka.gbk -o test_predictions
```

2. Compare the the 3Di sequences to the structure database with Foldseek using `phold compare`. This does not utilise a GPU. 

```
phold comapre --predictions_dir test_predictions -o test_output_phold -t 8 -d phold_structure_foldseek_db/
```

# Output

* The primary outputs are `final_cds_predictions.tsv`, where you have detailed annotation information on every CDS and `phold.gbk`, which contains a GenBank format file including these annotations, and any other genomic features (tRNA, CRISPR repeats, tmRNAs) from the `pharokka` input file.


## Usage

```
Usage: phold [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  citation  Print the citation(s) for this tool
  compare   Runs phold compare (Foldseek)
  createdb  Creates Foldseek DB from AA FASTA and 3Di FASTA input files
  predict   Runs phold predict (ProstT5)
  proteins  Runs phold proteins (ProstT5 on a multiFASTA input)
  remote    Runs phold predict using foldseek API and compare locally
  run       Runs phold predict (ProstT5) and comapare (Foldseek)
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
    * Probably will remove this



