#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by length variants"""

# Imports

from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys
import argparse
import os
import re
import warnings

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for splitting up the sequences in a fasta file according to a label in the sequence headers. Label format is assumed to be \";XXXX=label;\", where XXXX is a specified descriptor. Each unique value of label will generate a separate output file. \"label\" is assumed to be a string terminated by a semicolon, space or end-of-line. ")

parser.add_argument("input", help = "input file path", metavar = "FASTA")
parser.add_argument("-o","--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-d", "--descriptor",	help = "descriptor string used to identify the label", type = str, default = "sample", metavar = "XXXX")


# Class definitons

# Function definitions

if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.input))[0]
	
	# Make the output directory
	
	if not os.path.exists(args.output_directory):
		os.makedirs(args.output_directory)
	
	# Check for bad options
	
	# Set up parser
	
	regex = ";"+args.descriptor+"=([^\s;:]*)[;\s$]"
	
	# Set up dict for file handles
	
	files = {}
	
	# Output depending on options
	
	seqn = 0
	filen = 0
	
	with open(args.input) as infasta:
		
		for head, seq in SimpleFastaParser(infasta):
			seqn = seqn + 1
			m = re.search(regex, head)
			sample = 'unknown'
			if m:
				sample = m.group(1)
			else:
				warnings.warn("No label can be parsed from header " + " head")
			
			if not sample in files:
				filen = filen+1
				files[sample] = open(args.output_directory + "/" + sample + ".fasta", "w")
			
			files[sample].write(">%s\n%s\n" % (head, seq))
			
			sys.stdout.write("\r%i sequences split into %i files" % (seqn,filen))
			sys.stdout.flush()
			
		
	sys.stdout.write("\n")
	
	# Close files
	for name, handle in files.items():
		handle.close()
	


