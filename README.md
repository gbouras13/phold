# phold - phage annotation using protein structural homology

`phold` is sensititve annotation tool for bacteriophage genomes and metagenomes using protein strucutal homology. 

`phold` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to translate protein amino acid sequences to the 3Di token alphabet used by [foldseek](https://github.com/steineggerlab/foldseek). Foldseek is then used to search these against a database of 803k phage protein structures mostly predicted using [Colabfold](https://github.com/sokrypton/ColabFold). 

Alternatively, you can specify protein structures that you have pre-computed for your phage(s) instead of using ProstT5.

# Table of Contents

- [phold - phage annotation using protein structural homology](#phold---phage-annotation-using-protein-structural-homology)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Output](#output)
  - [Usage](#usage)
- [Citation](#citation)

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
  citation          Print the citation(s) for this tool
  compare           Runs Foldseek vs phold db
  createdb          Creates foldseek DB from AA FASTA and 3Di FASTA input...
  predict           Uses ProstT5 to predict 3Di tokens - GPU recommended
  proteins-compare  Runs Foldseek vs phold db on proteins input
  proteins-predict  Runs ProstT5 on a multiFASTA input - GPU recommended
  remote            Uses foldseek API to run ProstT5 then foldseek locally
  run               phold predict then comapare all in one - GPU recommended
```

* On large datasets, it is best to run `phold predict` (with GPU) followed by `phold compare` (CPU)


# Citation

`phold` is a work in progress, a preprint will be coming hopefully soon - if you use it please cite the github repository https://github.com/gbouras13/phold for now.

Please be sure to cite the following core dependencies and PHROGs database:

* [Foldseek](https://github.com/steineggerlab/foldseek) [van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist C, Söding J, and Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)](https://www.nature.com/articles/s41587-023-01773-0)
* [ProstT5](https://github.com/mheinzinger/ProstT5) [Michael Heinzinger, Konstantin Weissenow, Joaquin Gomez Sanchez, Adrian Henkel, Martin Steinegger, Burkhard Rost. ProstT5: Bilingual Language Model for Protein Sequence and Structure. bioRxiv doi:10.1101/2023.07.23.550085 (2023)](https://www.biorxiv.org/content/10.1101/2023.07.23.550085v1)
* [Colabfold](https://github.com/sokrypton/ColabFold) [Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold: Making protein folding accessible to all. Nature Methods (2022) doi: 10.1038/s41592-022-01488-1 ](https://www.nature.com/articles/s41592-022-01488-1)
* [PHROGs](https://phrogs.lmge.uca.fr) [Terzian P., Olo Ndela E., Galiez C., Lossouarn J., Pérez Bucio R.E., Mom R., Toussaint A., Petit M.A., Enault F., "PHROG : families of prokaryotic virus proteins clustered using remote homology", NAR Genomics and Bioinformatics, (2021) https://doi.org/10.1093/nargab/lqab067](https://doi.org/10.1093/nargab/lqab067)

Please also consider citing these supplementary databases where relevant:

* [CARD](https://phrogs.lmge.uca.fr) [Alcock B.P. et al, CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database Nucleic Acids Research (2022) https://doi.org/10.1093/nar/gkac920](https://doi.org/10.1093/nar/gkac920)
* [VFDB](http://www.mgc.ac.cn/VFs/main.htm) [Chen L., Yang J., Yao Z., Sun L., Shen Y., Jin Q., "VFDB: a reference database for bacterial virulence factors", Nucleic Acids Research (2005) https://doi.org/10.1093/nar/gki008](https://doi.org/10.1093/nar/gki008)
* [Defensefinder](https://defensefinder.mdmlab.fr) [ F. Tesson,  R. Planel, A. Egorov, H. Georjon,  H. Vaysset,  B. Brancotte,  B. Néron,  E. Mordret,  A Bernheim,  G. Atkinson,  J. Cury. A Comprehensive Resource for Exploring Antiphage Defense: DefenseFinder Webservice, Wiki and Databases. bioRxiv (2024) https://doi.org/10.1101/2024.01.25.577194](https://doi.org/10.1101/2024.01.25.577194)
* [acrDB](https://bcb.unl.edu/AcrDB/) - please cite the original acrDB database paper [Le Huang, Bowen Yang, Haidong Yi, Amina Asif, Jiawei Wang, Trevor Lithgow, Han Zhang, Fayyaz ul Amir Afsar Minhas, Yanbin Yin, AcrDB: a database of anti-CRISPR operons in prokaryotes and viruses. Nucleic Acids Research (2021) https://doi.org/10.1093/nar/gkaa857](https://doi.org/10.1093/nar/gkaa857) AND the paper that generated the structures for these protein used by `phold` [Harutyun Sahakyan, Kira S. Makarova, and Eugene V. Koonin. Search for Origins of Anti-CRISPR Proteins by Structure Comparison. The CRISPR Journal (2023)](https://doi.org/10.1089/crispr.2023.0011)


