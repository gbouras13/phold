# History

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