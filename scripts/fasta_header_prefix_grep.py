#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
    usage = 'python3 fasta_header_prefix_grep.py ...'
    parser = argparse.ArgumentParser(description='script to extract all records beginning with a string', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
    parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
    parser.add_argument('-p', '--prefix', action="store", help='string the header must contain',  required=True)

    args = parser.parse_args()

    return args

args = get_input()

# Create an output file to write the matching records
with open(args.outfile, "w") as output_handle:
    # Iterate through the input multi-FASTA file
    for record in SeqIO.parse(args.infile, "fasta"):
        # Check if the record's ID starts with the specified prefix
        if record.id.startswith(args.prefix):
            # Write the matching record to the output file
            SeqIO.write(record, output_handle, "fasta")
