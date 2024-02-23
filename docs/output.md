# `phold` Output Files

## Main Outputs

* `_aa.fasta` which will hold all amino acid sequences of predicted CDSs
* `_3di.fasta` which will hold all Foldseek 3Di sequences of predicted CDSs as predicted by ProstT5
* `_.gbk` which will contain a Genbank format file of your phage(s) with all annotations
* `_all_cds_functions.tsv` which includes for each contig:
  *  total CDS counts
  *  total CDS counts of each PHROG category 
  *  total CDS counts for CARD, VFDB, Defensefinder and ACR databases
*  `_per_cds_predictions.tsv` which contains detailed information for each CDS

Columns of this output include:

*  `phrog` which gives the tophit PHROG group
*  `function` which gives the PHROG category associated with the top hit
*  `product` which gives the detailed protein annotation
*  `bitscore` until `tLen` - these are taken from the foldseek output for the top hit protein
*  `function_with_highest_bitscore_proportion`  and all columns with `bitscore_proportion` are related to probabilistic annotation
   *  This is a work in progress
   *  For each query CDS, these columns give the proportion of the total Foldseek bitscore corresponding to that PHROG category of all hits for that CDS 
   *  For example, if a CDS had 10 Foldseek hits all with bitscore 100, and 9 of these hits were to 'tail' and 1 was to lysis, the value for 'tail_bitscore_proportion' would be 0.9 and the value for 'lysis_bitscore_proportion' would be 0.1.
   *  The intention of these columns are to indicate CDS that have more uncertainty regarding their annotation - whether they possibly have multiple domains, or simply because they are homologous to multiple different PHROGs - so you should probably have a closer look at them


## Supplementary Outputs

* `sub_db_tophits` which contains:
  * `acr_cds_predictions.tsv` - contains all CDS with top hit to the [Anti-CRISPR database](https://bcb.unl.edu/AcrDB/)
  * `defensefinder_cds_predictions.tsv` - contains all CDS with top hit to the [Defensefinder database](https://defensefinder.mdmlab.fr)
  * `card_cds_predictions.tsv` - contains all CDS with top hits to the [CARD AMR database](https://card.mcmaster.ca)
  * `vfdb_cds_predictions.tsv` - contains all CDS with top hits to the [VFDB database](vfdb_cds_predictions.tsv)
*  `_prostT5_3di_mean_probabilities.csv` - contains the mean ProstT5 probability score for each CDS. These are equivalent to the probability of how similar the overall ProstT5 3Di sequence is predicted to be compared to its Alphafold2 baseline
*  `_prostT5_3di_all_probabilities.json` - contains the ProstT5 probabilities for each residue for each CDS, in the json format
