#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Back-translate alignment of amino acids in a fasta to aligned nucleotide sequences using unaligned nucleotide sequences as a guide"""

# Imports
from Bio import SeqIO
import argparse
import sys

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for converting aligned amino acid sequences and the source unaligned nucleotide sequences into the equivalent aligned nucleotide sequences. All sequences are back-translated in the forward direction, using the same reading frame and translation table. Amino acid alignment should be supplied on STDIN, unaligned nucleotides should be supplied as a file path. Unaligned and aligned entries should be in exactly the same order. Results are written to STDOUT")

parser.add_argument("nucpath", help = "path to corresponding unaligned nucleotides", metavar = "NUCPATH")
parser.add_argument("table", help = "translation table number, required", choices = range(1,33), type = int, metavar = "TABLE")

parser.add_argument("-i","--input", help = "amino acid alignment", type = str)
parser.add_argument("-r","--reading_frame", help = "reading frame, default frame 1", type = int, choices = [1,2,3], default = 1)

# Function definitions

# Main

if __name__ == "__main__":
	
	# Get options
	
	#args = parser.parse_args(['-r', '1', '-i', 'aln_aa.fa', 'test.fa', '5'])
	args = parser.parse_args()
	
	rf0 = args.reading_frame-1
	
	# Read amino acid alignment
	if(args.input):
		aaa_records = SeqIO.parse(args.input, "fasta")
	else:
		aaa_records = SeqIO.parse(sys.stdin, "fasta")
	
	# Read nucleotides
	nuc_records = SeqIO.parse(args.nucpath, "fasta")
	
	# Compare lengths of generators
	#sys.exit("Error: number of sequences in amino acid alignment (", len(aaa_records), ") does not match number of sequences in unaligned nucleotides file (", len(nuc_records), ")")
	
	# Find gap positions and convert nucleotide data
	n = 0
	nuc_aln = list()
	for aaa_rec, nuc_rec in zip(aaa_records, nuc_records):
		#nuc_rec, aaa_rec = list(zip(aaa_records, nuc_records))[0]
		
		sys.stderr.write("Comparing AA %s %s with NT %s %s\n" % (aaa_rec.id, aaa_rec.seq, nuc_rec.id, nuc_rec.seq))
		
		n += 1
		
		# Run checks
		if(aaa_rec.id != nuc_rec.id):
			sys.exit("Error: sequence identifiers do not match for sequence " + str(n) + " (alignment: " + aaa_rec.id + " nucleotides: " + nuc_rec.id)
		
		# Add gaps for incomplete ends
		remainder = len(nuc_rec.seq[(args.reading_frame-1):]) % 3
		nuc_rec.seq = nuc_rec.seq + '-'*(3-remainder) if(remainder != 0) else nuc_rec.seq
		
		# Find gap positions
		aa_gaps = [i for i, b in enumerate(aaa_rec.seq) if b == "-"]
		
		# Convert to nucleotide positions, offsetting by reading frame
		#nuc_gaps = [b+rf0 for i in aa_gaps for b in range(i*3, i*3+2)]
		nuc_gaps = [a*3+rf0 for a in aa_gaps]
		
		# Insert gaps
		for g in nuc_gaps:
			nuc_rec.seq = nuc_rec.seq[:g] + "---" + nuc_rec.seq[g:]
		
		lendiff = len(str(aaa_rec.seq)) * 3 - (len(nuc_rec.seq) + rf0)
		if( lendiff != 0 ):
			if( lendiff < 3 ):
				nuc_rec.seq = nuc_rec.seq + ("-" * lendiff)
			else:
				sys.exit("Error: number of bases / amino acids do not correspond in sequence number " + str(n) + " " + aaa_rec.id + " (alignment)")
		
		nuc_aln.append(nuc_rec)
	
	# Write nucleotides
	SeqIO.write(nuc_aln, sys.stdout, "fasta")
