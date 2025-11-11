<p align="center">
  <img src="phold_logo.png" alt="phold Logo" height=200>
</p>

`phold` is a sensitive annotation tool for bacteriophage genomes and metagenomes using protein structural homology. 

`phold` uses the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model to rapidly translate protein amino acid sequences to the 3Di token alphabet used by [Foldseek](https://github.com/steineggerlab/foldseek). Foldseek is then used to search these against a database of over 1.36 million phage protein structures mostly predicted using [Colabfold](https://github.com/sokrypton/ColabFold). 

Alternatively, you can specify protein structures that you have pre-computed for your phage(s) instead of using ProstT5 with `phold compare`.

The `phold` databse consists of over 1.36 million protein structures generated using [Colabfold](https://github.com/sokrypton/ColabFold) and [ESMFold](https://github.com/facebookresearch/esm) from the following databases:

* [PHROGs](https://phrogs.lmge.uca.fr) — 441,177 de-deuplicated proteins. Proteins over 3000AA were chunked into equal components such that each fragment was under 3k.
* [enVhogs](http://envhog.u-ga.fr/envhog/) — 562,369 proteins from the ~2.2M ENVHOGs <3000AA that were assigned a PHROG.
* [efam](https://doi.org/10.1093/bioinformatics/btab451) - 233,181 efam "extra conservative" proteins that were assigned a PHROG.
* [DGRs](https://doi.org/10.1038/s41467-021-23402-7) - 12,683 extra diversity generating element proteins from [Roux et al 2021](https://doi.org/10.1038/s41467-021-23402-7).
* [VFDB](http://www.mgc.ac.cn/VFs/main.htm) —  27,823 structures of putative bacterial virulence factors from the VFDB.
* [CARD](https://card.mcmaster.ca) —  4,804 structures of anitbiotic resistant proteins from CARD.
* [acrDB](https://bcb.unl.edu/AcrDB/) -  3,652 anti-crispr proteins predicted in this [study](https://doi.org/10.1089/crispr.2023.0011).
* [Defensefinder](https://defensefinder.mdmlab.fr) - 408 monomer prokaryotic defense proteins.
* [Netflax](https://doi.org/10.1073/pnas.2305393120) - 7,152 toxin-antitoxin proteins.

# Google Colab Notebook

If you don't want to install `phold` locally, you can run it without any code using one [this Google Colab notebook](https://colab.research.google.com/github/gbouras13/phold/blob/main/run_pharokka_and_phold_and_phynteny.ipynb). 

Pharokka, Phold and Phynteny are complimentary tools and when used together, they substantially increase the annotation rate of your phage genome. The below plot shows the annotation rate of different tools across 4 benchmarked datasets ((a) INPHARED 1419, (b) Cook, (c) Crass and (d) Tara - see the [Phold preprint]((https://www.biorxiv.org/content/10.1101/2025.08.05.668817v1)) for more information)

Specifically, the final Phynteny plots combine the benefits of annotation with Pharokka (with HMM, the second violin) followed by Phold (with structures, the fourth violin) followed by Phynteny

<p align="center">
  <img src="Pharokka_Phold_Phynteny.png" alt="pharokka plus phold plus phynteny" height=300>
</p>