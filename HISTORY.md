# History

1.2.2 (2026-01-31)
------------------
* Minor bugfix to make phold is compatible with transformers v5 and pandas v3

1.2.1 (2026-01-15)
------------------
* Minor bugfix to make sure the ordering of ProstT5 confidence outputs matches the input .faa (was length-sorted, introduced by batching in v1.2.0)


1.2.0 (2026-01-08)
------------------
* Improved ProstT5 3Di prediction throughput for  phold run, phold predict and phold proteins-predict due to smarter batching implmentations
* Addition of phold autotune subcommand to detect an appropriate --batch_size for your hardware
* You can also use --autotune with phold run, phold predict and phold proteins-predict to automatically detect and use the optimal --batch_size (only recommended for large datasets with thousands of proteins)
* Manuscript now published: https://doi.org/10.1093/nar/gkaf1448

1.1.0 (2025-11-07)
------------------

* Integration with suvtk to make to it easier to submit Pharokka and Phold annotated genomes to Genbank - thanks to @LanderDC for suvtk and integration. See https://github.com/gbouras13/phold?tab=readme-ov-file#genbank-submission for more details
* Adds --restart parameter to complete large phold compare jobs #79


1.0.1 (2025-08-12)
------------------

* Minor bug fixes including
* #98 adding backup DB downloading if aria2c is not available https://github.com/gbouras13/phold/issues/98
* Adds check to ensure no colons (":") are the headers for FASTA/GenBank inputs (Phold uses delimiter) #86 . Please remove all colons from headers.
* Adds logic to handle GenBank formatted headers that Foldseek automatically trims to only the accession https://github.com/steineggerlab/foldseek/blob/cfb431e98abcc5bcc49950285211d3723b47dc94/lib/mmseqs/src/commons/Util.cpp#L138 #86


1.0.0 (2025-08-06)
------------------

Major Phold release to go with the preprint. For more details, see the preprint and updated documentation.

You will need to re-install the updated Phold search database with phold install to be compatible with v1.0.0

Major Changes

* Phold search database has been modified, filtered and curated to contain 1,363,704 proteins structures with functional labels (see https://zenodo.org/records/16741548). In particular, since the previous release of Phold, the enVhogs were re-clustered and re-labelled by the authors of that work. This release contains the updated enVhog structures.
* We additionally make available a larger database containing 3,166,602 structures (i.e. the Phold search database plus an extra 1.8M efam and enVhog proteins without PHROG assignment or functional label) to download using phold install --extended_db. Using this database provides marginally fewer functional annotations and takes longer than using the default Phold search database, so is not recommended for functional annotation, but finds more hits (i.e. including to unknown function proteins) overall, so may be of interest for viral identification tasks.
PHROG functional labels have been updated for 2,798 PHROGs using manual curation informed by structural similarity searches. See the preprint for more details. The updated annotations are available in the phold database under phold_annots.tsv
* Phold search database is no longer pre-clustered, as it was shown not to significantly differ in terms of sensitivity and runtime from unclustered for the updated database.
* Phold supports Foldseek-GPU acceleration for NVIDIA GPUs using `--foldseek_gpu`. Note that it is still ideal to run Phold with multiple CPU-threads (e.g. -t 8 or however many threads you have available), as GPU acceleration only accelerates and improves the prefilter of Foldseek.
* Phold supports custom user-specified Foldseek databases with `--custom_db`.
* Phold adds high, medium and low confidence annotation heuristics to guide the user (especially users from wet-lab backgrounds or without much understanding of protein structural alignment metrics) as to what annotations they should trust with a very high degree of confidence, and which they should prioritise for manual curation. See the documentation for more.
* Phold will now mask all residues below 25 ProstT5 Confidence by default (can be varied with `--mask_threshold`), as this was shown to increase annotation performance compared to no masking.
* If you only want to annotate hypothetical proteins from Pharokka to save runtime and resource usage, you can use `--hyps`
* You can run Phold with fine-tuned ProstT5 models using `--finetune` (phage finetuned ProstT5 encoder and phage fine-tuned CNN) or `--vanilla` (phage finetuned ProstT5 encoder and vanilla PDB-based CNN). Annotation performance with these do not dramatically differ with the default ProstT5 (see the preprint), but may be of interest to some users of Phold.
* Fix bug with `phold predict` with `--cpu`, where ProstT5 would use all threads by default https://github.com/gbouras13/phold/issues/60 Thanks @valentynbez 
* Fix bug with `phold compare` with `--structures`. If there were additional structures in the `--structure_dir` not found in the input CDS and `--filter_structures` was not specified, phold would crash if there was a Foldseek hit to the extra structures. Thanks Nikolas Basler
* Fixes bug that would occur if the ID/locus tag of the features in the input genbank were longer than 54 characters - Phold would complete but introduce a space into CDSids and FASTA headers (#67) 

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