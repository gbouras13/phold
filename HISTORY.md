# History

0.2.1 (2024-07-13)
------------------

* Fix bug with `phold predict` with `--cpu`, where ProstT5 would use all threads by default https://github.com/gbouras13/phold/issues/60 @valentynbez 
* Fix bug with `phold compare` with `--structures`. If there were additional structures in the `--structure_dir` not found in the input CDS and `--filter_structures` was not specified, phold would crash if there was a Foldseek hit to the extra structures.

0.2.0 (2024-07-13)
------------------

**You will need to re-install the updated phold database for v0.2.0 using `phold install`**
**You will also need to upgrade Foldseek to v9.427df8a**

v0.2.0 is a large update adding:

* Improved sensitivity and faster runtime for the `foldseek` search. This is achieved by clustering the Phold database at `--min-seq-id 0.3 -c 0.8` and creating a cluster db before running with `foldseek` which significantly improves runtime
    * Overall, just over 1.1M structures are clustered into around 372k clusters 
* `--cluster-search 1` parameter is added to `foldseek search` to search against the cluster representatives first and then within each cluster, which increases sensitivity and reduces resource usage compared to `phold v0.1.4`
* Changed default `--max_seqs` from 1000 to 10000 to improve sensitivity at little resource usage cost
* Phold database is expanded adding:
    * Extremely conservative high confidence [efam](https://doi.org/10.1093/bioinformatics/btab451) proteins with hits to PHROGs.
    * 95% dereplicated diversity-generating retroelements (DGRs) from [Roux et al](https://www.nature.com/articles/s41467-021-23402-7).
    * 7153 netflax toxin-antitoxin system proteins from [Ernits et al](https://doi.org/10.1073/pnas.2305393120).
* Adds `--ultra_sensitive` flag which turns off Foldseek prefiltering for maximum sensitivity. Recommended for small datasets/single phages only.
    * This passes the `--exhaustive-search` parameter to `foldseek search`
* Adds the ability to save ProstT5 embeddings with `--save_per_residue_embeddings` and `--save_per_protein_embeddings`
* Adds `.cif` support (e.g. from Alphafold3 server) for structures, not just `.pdb` file format
* Removes some experimental parameters from v0.1.4 (`--split` etc)

Breaking CLI parameter changes

* `--pdb` has changed to `--structures`
* `--pdb_dir` has changed to `--structure_dir`
* `--filter_pdbs` has changed to `--filter_structures`

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