# `phold` Output Files

## Main Outputs

* `_aa.fasta` which will hold all amino acid sequences of predicted CDSs
* `_3di.fasta` which will hold all Foldseek 3Di sequences of predicted CDSs as predicted by ProstT5
* `_.gbk` which will contain a Genbank format file of your phage(s) with all annotations
* `_all_cds_functions.tsv` which includes for each contig:
    * Total CDS counts
    * Total CDS counts of each PHROG category 
    * Total CDS counts for CARD, VFDB, Defensefinder and ACR databases
*  `_per_cds_predictions.tsv` which contains detailed information for each CDS

The columns of `_per_cds_predictions.tsv` include:

*  `phrog` which gives the tophit PHROG group
*  `function` which gives the PHROG category associated with the top hit
*  `product` which gives the detailed protein annotation
*  `bitscore` until `tLen` - these are taken from the foldseek output for the top hit protein
*  `function_with_highest_bitscore_proportion`  and all columns with `bitscore_proportion` are related to probabilistic annotation
    * For each query CDS, these columns give the proportion of the total Foldseek bitscore corresponding to that PHROG category of all hits for that CDS 
    * For example, if a CDS had 10 Foldseek hits all with bitscore 100, and 9 of these hits were to 'tail' and 1 was to lysis, the value for 'tail_bitscore_proportion' would be 0.9 and the value for 'lysis_bitscore_proportion' would be 0.1.
    * The intention of these columns are to indicate CDS that have more uncertainty regarding their annotation - whether they possibly have multiple domains, or simply because they are homologous to multiple different PHROGs - so you should probably have a closer look at them
* `annotation_confidence` - All annotations are denoted as high, medium or low confidence on the basis of the following heuristics to make the output of Phold more interpretable for users, particularly those without understanding of protein structural alignment methods:
    * High confidence hits are where both the query and target proteins have at least 80% reciprocal alignment coverage, along with either (i) greater than 30% amino acid sequence identity (suggesting the hit is in the light zone of sequence homology27) or (ii) query mean ProstT5 confidence of at least 60% (suggesting a very good quality ProstT5 3Di prediction) or (iii) an alignment E-value < 1e-10. 
    * Medium confidence hits are where either the query or target protein hit has at least 80% coverage along with either (i) greater than 30% amino acid sequence identity or (ii) ProstT5 confidence between 45% and 60% (suggesting a good quality ProstT5 3Di prediction) and (iii) an E-value < 1-e05. 
    * Low confidence hits are all other hits below the specified E-value that do not fit those thresholds (i.e. hits with low coverage, low amino acid sequence identity, low ProstT5 confidence, or hits near the E-value threshold of 0.001). 
    * If Phold is run with user provided input structures instead of ProstT5, the heuristics are identical other than the ProstT5 confidence criteria. 
    * In this case, Phold will also output alignment template modeling score (TM-score) and local distance difference test (LDDT) values from Foldseek, which may also guide the user in assessing annotation quality. 
    * To be clear, 'low' quality annotations do not mean they are wrong, only that these annotations should not be trusted as much as high (extremely trustworthy) or medium (very trustworthy) annotations. They probably need a bit more manual annotation than the others - something that any good phage biologist should do for phage annotations anyway given how hard phages are to annotate!

## Supplementary Outputs

* `sub_db_tophits` which contains:
    * `acr_cds_predictions.tsv` - contains all CDS with top hit to the [Anti-CRISPR database](https://bcb.unl.edu/AcrDB/)
    * `defensefinder_cds_predictions.tsv` - contains all CDS with top hit to the [Defensefinder database](https://defensefinder.mdmlab.fr)
    * `card_cds_predictions.tsv` - contains all CDS with top hits to the [CARD AMR database](https://card.mcmaster.ca)
    * `vfdb_cds_predictions.tsv` - contains all CDS with top hits to the [VFDB database](vfdb_cds_predictions.tsv)
    * `netflax_cds_predictions.tsv` - contains all CDS with top hits to NetFlax [toxin-antitoxins](https://doi.org/10.1073/pnas.2305393120).
*  `_prostT5_3di_mean_probabilities.csv` - contains the mean ProstT5 probability score for each CDS. These are equivalent to the probability of how similar the overall ProstT5 3Di sequence is predicted to be compared to its Alphafold2 baseline
*  `_prostT5_3di_all_probabilities.json` - contains the ProstT5 probabilities for each residue for each CDS, in the json format

## Optional Outputs
*  `_embeddings_per_protein.h5` - contains the ProstT5 embeddings for each protein in the h5 format
*  `_embeddings_per_residue.h5` - contains the ProstT5 embeddings for each residue in the h5 format
*  `phold_custom_database_hits.tsv` - contains Phold custom database hits is `--custom_db` is specified
