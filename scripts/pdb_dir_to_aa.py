#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import PDB
import os

def get_input():
    usage = 'python3 pdb_dir_to_aa.py ...'
    parser = argparse.ArgumentParser(description='script to extract AAs from a directory of PDB files and write as a multiFASTA', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--indir', action="store", help='Input directory of PDBs',  required=True)
    parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
    args = parser.parse_args()

    return args



def extract_amino_acid_sequence_from_file(pdb_file):

    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure('protein', pdb_file)

    # Initialize an empty sequence
    sequence = ''

    # Iterate over all models in the structure (usually there is only one)
    for model in structure:
        # Iterate over all chains in the model
        for chain in model:
            # Iterate over all residues in the chain
            for residue in chain:
                # Check if the residue is an amino acid
                if PDB.is_aa(residue):
                    # Get the three-letter code of the amino acid
                    aa_code = PDB.Polypeptide.protein_letters_3to1(residue.get_resname())
                    
                    # Append the amino acid code to the sequence
                    sequence += aa_code

    return sequence


def main():

    args = get_input()

    # seq dic
    pdb_sequences = {}

    # List pdb all files in the directory
    pdb_files = [f for f in os.listdir(args.indir) if f.endswith('.pdb')]

    for pdb_file in pdb_files:
        pdb_file_path = os.path.join(args.indir, pdb_file)
        amino_acid_sequence = extract_amino_acid_sequence_from_file(pdb_file_path)
        phrog_name = pdb_file.replace(".pdb", "")
        pdb_sequences[phrog_name] = amino_acid_sequence


    ## write the CDS to file
    fasta_aa = args.outfile
    with open(fasta_aa, "w+") as out_f:
        for phrog_id, aa_seq in pdb_sequences.items():
            out_f.write(f">{phrog_id}\n")
            out_f.write(f"{pdb_sequences[phrog_id]}\n")


if __name__ == "__main__":
    main()


