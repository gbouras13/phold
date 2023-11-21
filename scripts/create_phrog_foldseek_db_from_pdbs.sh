#!/bin/bash

# chopped ones

find phrog_prediction_chopped/ -name '*unrelaxed_rank_001*' -exec cp {} test_pdb/ \;

# for the chopped ones
cd test_pdb
find . -name 'phrog*_000.pdb' -maxdepth 1 -print0 | xargs -0 -n 10 rename 's/phrog_(\d+)_(\d+)_unrelaxed_rank_001_alphafold2_ptm_model_\d+_seed_000.pdb/phrog_$1_$2.pdb/'

cd test_pdb ..

## for all the non chopped ones
find phrog_prediction/ -name '*unrelaxed_rank_001*' -exec cp {} test_pdb/ \;
cd test_pdb

# renames the extra stuff
find . -name 'phrog*.pdb' -maxdepth 1 -print0 | xargs -0 -n 10 rename 's/phrog_(\d+)_unrelaxed_rank_001_alphafold2_ptm_model_\d+_seed_000.pdb/phrog_$1.pdb/'
cd test_pdb ..

# to actually create the db
mkdir PHROG_Foldseek_PDB_db
foldseek createdb test_pdb/ PHROG_Foldseek_PDB_db/all_phrogs_pdb
