#!/bin/bash

# generally speaking

<input.fasta>

# phanotate

pharokka.py -i <input.fasta > -o <output_dir> -d <db> -t 16 

# prodigal

pharokka.py -i <input.fasta > -o <output_dir> -d <db> -t 16 -g prodigal


# run colabfold. Whether locally or not needs to be on the phanotate.faa/prodigal.faa file
# localcolabfold (for me)
colabfold_batch --amber --num-relax 1 --num-recycler 3 --relax-max-iterations 5 --use-gpu-relax <phanotate.faa> <colabfold_output_dir>

# phold colabfold

phold predict

phold compare 
