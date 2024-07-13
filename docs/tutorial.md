# `phold` Tutorial

* This tutorial assumes you have [mamba](https://github.com/conda-forge/miniforge) installed and the correct channels available. Please see the [install page](https://phold.readthedocs.io/en/latest/install/) for more details. You could also use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* This tutorial uses _Stenotrophomonas_ phage SMA6, accession NC_043029.

## Step 1 Get the `phold` repository and test data

```bash
git clone "https://github.com/gbouras13/phold.git"
cd phold
```

* The test genome `NC_043029.fasta` should be in the `tests/test_data`

## Step 2 Installing and Running `Pharokka`

* Feel free to skip this step if you only want phage CDS annotations from `phold` - Pharokka also gives you tRNA, tmRNAs and CRISPRs and other summary features `phold` lacks (at least for now).


* To install and run pharokka (change `-t 8` to the number of available threads)

```bash
mamba create -n pharokkaENV pharokka
conda activate pharokkaENV
install_databases.py -o pharokka_db
pharokka.py -i tests/test_data/NC_043029.fasta -d pharokka_db -o NC_043029_pharokka_output -t 8 --fast
conda deactivate
```

## Step 3 Installing `phold`

* To install `phold` with mamba from bioconda (assuming you have an NVIDIA GPU available):

```bash
mamba create -n pholdENV -c conda-forge -c bioconda phold pytorch=*=cuda*
conda activate pholdENV
phold install
```

* If you do not have a GPU available, remove `pytorch=*=cuda*`
* For more installation options, see the [installation documentation](https://phold.readthedocs.io/en/latest/install/).

## Step 4 Running `phold`

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029
```

* If you do not have a GPU available:

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029 --cpu
```

* If you skipped step 2, replace `NC_043029_pharokka_output/pharokka.gbk` with `tests/test_data/NC_043029.fasta`

## Step 5 Running `phold plot`

* `phold` can generate Circos plot of your phage(s)
* The plot will be saved in the `NC_043029_phold_plots` directory. See the [documentation](https://phold.readthedocs.io/en/latest/run/#phold-plot) for more parameter details
* `phold plot` provides .png and .svg outputs

```bash
phold plot -i NC_043029_phold_output/NC_043029.gbk -o NC_043029_phold_plot -t '${Stenotrophomonas}$ Phage SMA6'
```

![Image](NC_043029.png)



