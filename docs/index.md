<p align="center">
  <img src="phold_logo.png" alt="phold Logo" height=200>
</p>

`phold` is a sensitive annotation tool for bacteriophage genomes and metagenomes using protein structural homology. 

`phold` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to rapidly translate protein amino acid sequences to the 3Di token alphabet used by [Foldseek](https://github.com/steineggerlab/foldseek). Foldseek is then used to search these against a database of over 1 million phage protein structures mostly predicted using [Colabfold](https://github.com/sokrypton/ColabFold). 

Alternatively, you can specify protein structures that you have pre-computed for your phage(s) instead of using ProstT5 with `phold compare`.

The `phold` databse consists of over 1 million protein structures generated using [Colabfold](https://github.com/sokrypton/ColabFold) from the following databases:

* [PHROGs](https://phrogs.lmge.uca.fr) — 440549 de-deuplicated proteins. Proteins over 1000AA were chunked into 1000AA components.
* [ENVHOGs](http://envhog.u-ga.fr/envhog/) — 315k representative proteins from the 2.2M ENVHOGs that have a PHROG function that is not 'unknown function'. Proteins over 1000AA were chunked into 1000AA components.
* [EFAM](https://doi.org/10.1093/bioinformatics/btab451) - 262k efam "extra conservative" proteins. Proteins over 1000AA were chunked into 1000AA components.
* [DGRs](https://doi.org/10.1038/s41467-021-23402-7) - 12683 extra diversity generating element proteins from [Roux et al 2021](https://doi.org/10.1038/s41467-021-23402-7).
* [VFDB](http://www.mgc.ac.cn/VFs/main.htm) — over 28k structures of putative bacterial virulence factors from the VFDB.
* [CARD](https://card.mcmaster.ca) — nearly 5k structures of anitbiotic resistant proteins from CARD.
* [acrDB](https://bcb.unl.edu/AcrDB/) - nearly 3.7k anti-crispr proteins predicted in this [study](https://doi.org/10.1089/crispr.2023.0011).
* [Defensefinder](https://defensefinder.mdmlab.fr) - 455 monomer prokaryotic defense proteins.
* [Netflax](https://doi.org/10.1073/pnas.2305393120) - 7153 toxin-antitoxin proteins.

# Google Colab Notebook

If you don't want to install `phold` locally, you can run it without any code using one [this Google Colab notebook](https://colab.research.google.com/github/gbouras13/phold/blob/main/run_pharokka_and_phold_and_phynteny.ipynb)

