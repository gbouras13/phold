# History

0.2.0
------------------

Big update adding

* Improved sensitivity and faster runtime for the `foldseek` search. This is achieved by clustering the Phold database at `--min-seq-id 0.3` and creating a cluster db before running with `foldseek` which significantly improves runtime
    * **This means you will need to re-install the phold database updates for v0.2.0 using `phold install`**
    * Additionally, changed default `--max_seqs` from 1000 to 10000 to improve sensitivity
* Phold database is expanded adding:
    * Extremely conservative high confidence [efam](https://doi.org/10.1093/bioinformatics/btab451) proteins with hits to PHROGs 
    * 95% dereplicated diversity-generating retroelements (DGRs) from [Roux et al](https://www.nature.com/articles/s41467-021-23402-7)
* Adds `depolymerase` column in `_per_cds_predictions.tsv` output.
    * Every PHROG, ENVHOG and efam protein in the Phold database has been run through [DepoScope](https://github.com/dimiboeckaerts/DepoScope)  
* Adds `--ultra_sensitive` flag which turns off Foldseek prefiltering for maximum annotation sensitivity. Recommended for small datasets/single phages only.
    * This passes the `--exhaustive-search` parameter to `foldseek search`
* Adds the ability to save ProstT5 embeddings with `--save_per_residue_embeddings` and `--save_per_protein_embeddings`


0.1.4 (2024-03-26)
------------------

* Fixes #31 issue with older Pharokka genbank input (prior to v1.5.0) that lacked 'transl_table' field
    * All Pharokka genbank input prior to v1.5.0 will be transl_table 11 (it is before pyrodigal-gv was added)
* Fixes genbank parsing bug that would occur if the ID/locus tag of the features in the inout genbank were longer than 54 characters 

0.1.3 (2024-03-19)
------------------

* Adds compability with Apple Silicon (M1/M2/M3) GPUs
* Fixes memory issue for `phold plot` with many contigs

0.1.2 (2024-03-06)
------------------

* Fixes `phold compare` cds_id issue where input file was FASTA
* Fixes issues with `phold remote` where input file was FASTA
* Improved documentation with conda/mamba install

0.1.1 (2024-03-05)
------------------

* Restructuring for pip and conda installation
* No substantive changes to the code base

0.1.0 (2024-03-05)
------------------

* Initial beta release