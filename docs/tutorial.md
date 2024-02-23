# `phold` Tutorial

* This tutorial assumes you have [mamba](https://github.com/conda-forge/miniforge) installed and the correct channels available. Please see the [install page](https://phold.readthedocs.io/en/latest/install/) for more details.

## Step 1 Get the `phold` repository and test data

```bash
git clone "https://github.com/gbouras13/phold.git"
cd phold
```

* The test genome `NC_043029.fasta` should be in the `tests/test_data`

## Step 2 Installing and Running `Pharokka`

* Feel free to skip this step if you only want phage CDS annotations from `phold` - Pharokka also gives you tRNA, tmRNAs and CRISPRs, plots and other summary features `phold` lacks (at least for now).


* To install and run pharokka (change `-t 8` to the number of available threads)

```bash
mamba create -n pharokkaENV pharokka
conda activate pharokkaENV
install_databases.py -o pharokka_db
pharokka.py -i tests/test_data/NC_043029.fasta -d pharokka_db -o NC_043029_pharokka_output -t 8 --fast
conda deactivate
```

## Step 3 Installing `phold`

* To install phold from source, replace the pip step with `pip install -e .`

```bash
mamba create -n pholdENV foldseek pip
conda activate pholdENV
pip install phold 
phold install
```

## Step 4 Running `phold`

* If you skipped step 2, replace `NC_043029_pharokka_output/pharokka.gbk` with `tests/test_data/NC_043029.fasta`

* If you have a GPU available:

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029
conda deactivate
```

* If you do not have a GPU available:

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029 --cpu
conda deactivate
```