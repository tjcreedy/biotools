#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Translate nucleotide sequences in a fasta to amino acids"""

# Imports
from Bio import SeqIO
import argparse
import sys

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for translating nucleotide sequences in a multifasta, supplied on STDIN. All sequences are translated in the forward direction, using the same reading frame and translation table. Results are written to STDOUT")

parser.add_argument("table", help = "translation table number, required", choices = range(1,33), type = int, metavar = "TABLE")
parser.add_argument("-r","--reading_frame",help = "reading frame, default frame 1", type = int, choices = [1,2,3], default = 1)

# Function definitions


# Main

if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Read nucleotides
	nuc_records = SeqIO.parse(sys.stdin, "fasta")
	
	
	# Translate nucleotides
	aa_records = list()
	for nuc_rec in nuc_records:
		aa_rec = nuc_rec
		aa_rec.seq = nuc_rec.seq[(args.reading_frame-1):].translate(table = args.table)
		aa_records.append(aa_rec)
	
	# Write amino acids
	SeqIO.write(aa_records, sys.stdout, "fasta")
