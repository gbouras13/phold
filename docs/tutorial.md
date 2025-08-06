# `phold` Tutorial

* This tutorial assumes you have [conda](https://github.com/conda-forge/miniforge) installed and the correct channels available. Please see the [install page](https://phold.readthedocs.io/en/latest/install/) for more details. You could also use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
* This tutorial uses _Stenotrophomonas_ phage SMA6, accession NC_043029.

## Step 1 Get the `phold` repository and test data

```bash
git clone "https://github.com/gbouras13/phold.git"
cd phold
```

* The test genome `NC_043029.fasta` should be in the `tests/test_data`

## Step 2 Installing and Running `Pharokka`

* Feel free to skip this step if you only want phage CDS annotations from `phold` - Pharokka also gives you tRNA, tmRNAs and CRISPRs and other summary features `phold` lacks.

* To install and run pharokka (change `-t 8` to the number of available threads)

```bash
conda create -n pharokkaENV pharokka
conda activate pharokkaENV
install_databases.py -o pharokka_db
pharokka.py -i tests/test_data/NC_043029.fasta -d pharokka_db -o NC_043029_pharokka_output -t 8 --fast
conda deactivate
```

## Step 3 Installing `phold`

* To install `phold` with conda from bioconda (assuming you have an NVIDIA GPU available):

```bash
conda create -n pholdENV -c conda-forge -c bioconda phold pytorch=*=cuda*
conda activate pholdENV
phold install -t 8 --foldseek_gpu
```

* If you do not have an NVIDIA GPU available, remove `pytorch=*=cuda*` and `--foldseek_gpu`
* For more installation options including non-NVIDIA GPUs including Mac Mseries laptops, see the [installation documentation](https://phold.readthedocs.io/en/latest/install/).

## Step 4 Running `phold`

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029 --foldseek_gpu
```

* If you do not have a GPU available:

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029 --cpu
```

* If you have a non-NVIDIA GPU available (e.g. Mac Mseries) and have installed the correct compatible version of PyTorch:

```bash
phold run -i NC_043029_pharokka_output/pharokka.gbk -o NC_043029_phold_output -t 8 -p NC_043029
```

* If you skipped step 2, replace `NC_043029_pharokka_output/pharokka.gbk` with `tests/test_data/NC_043029.fasta`


## Step 4.5 Interpreting Phold Output

* The `annotation_confidence` column of is a good start to determine the quality of your phage annotations `_per_cds_predictions.tsv`. All annotations are either high, medium or low confidence based on the following heuristics:
    * High confidence hits are where both the query and target proteins have at least 80% reciprocal alignment coverage, along with either (i) greater than 30% amino acid sequence identity (suggesting the hit is in the light zone of sequence homology27) or (ii) query mean ProstT5 confidence of at least 60% (suggesting a very good quality ProstT5 3Di prediction) or (iii) an alignment E-value < 1e-10. 
    * Medium confidence hits are where either the query or target protein hit has at least 80% coverage along with either (i) greater than 30% amino acid sequence identity or (ii) ProstT5 confidence between 45% and 60% (suggesting a good quality ProstT5 3Di prediction) and (iii) an E-value < 1-e05. 
    * Low confidence hits are all other hits below the specified E-value that do not fit those thresholds (i.e. hits with low coverage, low amino acid sequence identity, low ProstT5 confidence, or hits near the E-value threshold of 0.001). 
    * If Phold is run with user provided input structures instead of ProstT5, the heuristics are identical other than the ProstT5 confidence criteria. 
    * In this case, Phold will also output alignment template modeling score (TM-score) and local distance difference test (LDDT) values from Foldseek, which may also guide the user in assessing annotation quality. 

To be clear, 'low' quality annotations do not mean they are wrong, only that these annotations should not be trusted as much as high (extremely trustworthy) or medium (very trustworthy) annotations. They probably need a bit more manual annotation than the others - something that any good phage biologist should do for phage annotations anyway given how hard phages are to annotate!

## Step 5 Running `phold plot`

* `phold` can generate Circos plot of your phage(s)
* The plot will be saved in the `NC_043029_phold_plots` directory. See the [documentation](https://phold.readthedocs.io/en/latest/run/#phold-plot) for more parameter details
* `phold plot` provides .png and .svg outputs

```bash
phold plot -i NC_043029_phold_output/NC_043029.gbk -o NC_043029_phold_plot -t '${Stenotrophomonas}$ Phage SMA6'
```

![Image](NC_043029.png)



