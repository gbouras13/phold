`phold` is sensititve annotation tool for bacteriophage genomes and metagenomes using protein strucutal homology. 

`phold` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to translate protein amino acid sequences to the 3Di token alphabet used by [foldseek](https://github.com/steineggerlab/foldseek). Foldseek is then used to search these against a database of 803k protein structures mostly predicted using [Colabfold](https://github.com/sokrypton/ColabFold). 

Alternatively, you can specify protein structures that you have pre-computed for your phage(s) instead of using ProstT5.

The `phold` databse consists of approximately 803k protein structures generated using [Colabfold](https://github.com/sokrypton/ColabFold) from the following databases:

* [PHROGs](https://phrogs.lmge.uca.fr) — 440549 de-deuplicated proteins. Proteins over 1000AA were chunked into 1000AA components.
* [ENVHOGs](http://envhog.u-ga.fr/envhog/) — 315k representative proteins from the 2.2M ENVHOGs that have a PHROG function that is not 'unknown function'. Proteins over 1000AA were chunked into 1000AA components.
* [VFDB](http://www.mgc.ac.cn/VFs/main.htm) — over 28k structures of putative bacterial virulence factors from the VFDB.
* [CARD](https://card.mcmaster.ca) — nearly 5k structures of anitbiotic resistant proteins from CARD.
* [acrDB](https://bcb.unl.edu/AcrDB/) - nearly 3.7k anti-crispr proteins predicted in this [study](https://doi.org/10.1089/crispr.2023.0011).
* [Defensefinder](https://defensefinder.mdmlab.fr) - 455 monomer prokaryotic defense proteins.


