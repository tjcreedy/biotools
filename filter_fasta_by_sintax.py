#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys
import argparse
import os
import re
import warnings

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for filtering a fasta file based on taxonomy in a SINTAX-formatted file. Supply a fasta and one or taxa, including the syntax level specifier (e.g. k:Eukaryota , o:Hymenoptera, etc). By default, sequences matching any of the input taxa will be output, and all others rejected. To select only sequences that match all of the filters, use the --all argument. To discard sequences matching the filters and retain the others, use the --reject argument.")

parser.add_argument("-f", "--fasta",		help = "fasta file path", metavar = "FASTA", type = str)
parser.add_argument("-s", "--sintax",		help = "sintax file path", metavar = "SINTAX", type = str)
parser.add_argument("-t", "--taxa",	help = "one or more taxa with level specifiers", metavar = 
"TAXON", nargs = '+')
parser.add_argument("-a", "--all",	help = "select sequences matching all filters (as opposed to any)", default = False, action = "store_true")
parser.add_argument("-r", "--reject",	help = "reject the selected taxa (as opposed to retaining)", default = False, action = "store_true")



# Class definitons

# Function definitions

if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Check for bad options
	
	# Parse filter
	
	selected = list()
	
	with open(args.sintax, "r") as insin:
		line = insin.readline()
		while line:
			line = line.strip()
			present = [f in line for f in args.taxa]
			if( (args.all and all(present)) or (not args.all and any(present)) ):
				selected.append(line.split("\t")[0])
			line = insin.readline()
	
	# Output depending on options
	
	with open(args.fasta, "r") as infasta:
		for head, seq in SimpleFastaParser(infasta):
			if head in selected:
				sys.stdout.write(">%s\n%s\n" % (head, seq) ) 
		
	
	


