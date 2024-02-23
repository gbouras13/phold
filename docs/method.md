# `phold` Method

## Step 1 ProstT5 3Di Inference

* `phold` begins by predicting the Foldseek 3Di tokens for every input protein using the [ProstT5](https://github.com/mheinzinger/ProstT5) protein language model
* Alternatively, this step is skipped if you specify pre-computed protein structures in the .pdb format using `--pdb_dir`

## Step 2 Foldseek Structural Comparison

* `phold` then creates a [Foldseek](https://github.com/steineggerlab/foldseek) database combining the AA and 3Di representations of each protein, and compared this to the `phold` database with Foldseek

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

* The `phold` database contains approximately 803 thousand phage protein structures generated mostly with Colabfold v1.5.3 from the following sources:
  * 766k [PHROG](https://phrogs.lmge.uca.fr) proteins
  * 3.6k [anti-CRISPR](https://bcb.unl.edu/AcrDB/) proteins
  * 455 [Defensefinder](https://defensefinder.mdmlab.fr) proteins
  * 28k [VFDB](http://www.mgc.ac.cn/VFs/main.htm) proteins (virulence factors)
  * 5k [CARD](https://card.mcmaster.ca) (AMR proteins)
