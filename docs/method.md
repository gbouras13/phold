# `phold` Method

<p align="center">
  <img src="img/phold_workflow.png" alt="phold method" height=300>
</p>

## Step 1 ProstT5 3Di Inference

* `phold` begins by predicting the Foldseek 3Di tokens for every input protein using the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model
* Alternatively, this step is skipped if you run `phold compare` and specify pre-computed protein structures in the .pdb or .cif formats using the `--structures` flag along with specifying the firectory containing the structures with `--structure_dir`

## Step 2 Foldseek Structural Comparison

* `phold` then creates a [Foldseek](https://github.com/steineggerlab/foldseek) database combining the AA and 3Di representations of each protein, and compares this to the `phold` database with Foldseek
* Alternatively, you can specify protein structures that you have pre-computed for your phage(s) instead of using ProstT5 with the parameter `--structures` and `--structure_dir`

## Step 3 Downstream Processing

* For each protein the following logic is conducted to select the top hit (if the input is not a Pharokka Genbank, the pharokka steps are skipped)

  * If there is at least 1 Foldseek hit that is non-hypothetical function, the top hit is considered to be that with the lowest evalue
  * Otherwise, if Pharokka non-hypothetical function hit exists in the Genbank input, take is considered the top hit
  * Otherwise, if only hypothetical proteins are found as foldseek hits, the hit with the lowest evalue (this PHROG's function may be discovered in the future)
  * Otherwise, if there are no foldseek hits and only a Pharokka hypothetical hit exists from the Genbank input, that is taken
  * Otherwise, if there are no foldseek hits and no Pharokka hit, then it is labelled ‘No_PHROG’ + hypothetical protein

## Step 4 Probabilistic Annotation

* **This is a work in progress**
* The aim is to reflect phage proteins where there exists uncertatinty of their function
* For each query CDS, all Foldseek hits above the Evalue threshold are consider
* The Foldseek bitscore sum is calculated for all 9 non-unknown PHROG functional categories
* This is divided by the total bitscore sum across all 9 categories
   *  For example, if a CDS had 10 Foldseek hits all with bitscore 100, and 9 of these hits were to 'tail' and 1 was to lysis, phold would consider the bitscore proportion value for tail to be 0.9 and the value for bitscore proportion value for lysis would be 0.1.
   *  The intention of these columns are to indicate CDS that have more uncertainty regarding their annotation - whether they possibly have multiple domains, or simply because they are homologous to multiple different PHROGs - so you should probably have a closer look at them
   *  This information can be found in `_per_cds_predictions.tsv`

## `phold` database

The `phold` databse consists of over 1 million protein structures generated mostly using [Colabfold](https://github.com/sokrypton/ColabFold) v1.5.3 (other than acrDB, Defensefinder and Netflax) from the following databases:

* [PHROGs](https://phrogs.lmge.uca.fr) — 440549 de-deuplicated proteins. Proteins over 1000AA were chunked into 1000AA components.
* [ENVHOGs](http://envhog.u-ga.fr/envhog/) — 315k representative proteins from the 2.2M ENVHOGs that have a PHROG function that is not 'unknown function'. Proteins over 1000AA were chunked into 1000AA components.
* [EFAM](https://doi.org/10.1093/bioinformatics/btab451) - 262k efam "extra conservative" proteins. Proteins over 1000AA were chunked into 1000AA components.
* [DGRs](https://doi.org/10.1038/s41467-021-23402-7) - 12683 extra diversity generating element proteins from [Roux et al 2021](https://doi.org/10.1038/s41467-021-23402-7).
* [VFDB](http://www.mgc.ac.cn/VFs/main.htm) — over 28k structures of putative bacterial virulence factors from the VFDB.
* [CARD](https://card.mcmaster.ca) — nearly 5k structures of anitbiotic resistant proteins from CARD.
* [acrDB](https://bcb.unl.edu/AcrDB/) - nearly 3.7k anti-crispr proteins predicted in this [study](https://doi.org/10.1089/crispr.2023.0011).
* [Defensefinder](https://defensefinder.mdmlab.fr) - 455 monomer prokaryotic defense proteins.
* [Netflax](https://doi.org/10.1073/pnas.2305393120) - 7153 toxin-antitoxin proteins.
