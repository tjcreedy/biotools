#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 09:52:39 2020

@author: thomas
"""

# Imports

import sys
import re
import argparse

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for correcting the format of headers in a genbank file to ensure proper parsing by BioPython. Supply a genbank file with one or more entries")

parser.add_argument("input", help = "input genbank file path", metavar = "file.gb")
parser.add_argument("-d","--default_division", help = "a three-letter code to be used as the default devision where the division code is missing", default = "UNK", metavar = "DIV")
parser.add_argument("-s", "--show_differences", help = "rather than the corrected genbank file, output old and new header lines", action = 'store_true')

# Class definitons

# Function definitions

def checkn(lis, name, ln, alt = None):
	n = len(lis)
	if n == 0: return alt
	if n == 1: return lis[0]
	if n > 1: sys.exit("Line %s has multiple matches for %s" % (str(ln), name))

# Main

if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	#args = parser.parse_args(['-d', 'INV', '-s', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-04-24_current/CCCP00043.gb'])
	
	# Open file and set up handler
	gb = open(args.input, 'r')
	
	ln = 0
	# Work through file
	for line in gb:
		#line = gb.readline()
		line = line.rstrip()
		ln += 1
		
		# Check if header line
		if re.match("^LOCUS", line):
			
			# Length
			length = re.findall('\d+ +(?:bp(?= ))', line)
			length = checkn(length, "sequence length", ln, '  bp')
			
			# Type
			seqtype = [ s for s in ['DNA', 'RNA', 'tRNA', 'mRNA', 'uRNA', 'snRNA', 'cDNA'] if ' ' + s + ' ' in line]
			seqtype =  '{:<7}'.format(checkn(seqtype, "sequence type", ln, ''))
			
			# Strand
			strand = [s for s in ['ss-', 'ds-', 'ms-'] if ' ' + s + ' ' in line]
			strand = checkn(strand, "sequence strand", ln, '   ')
			
			# Date
			date = re.findall('(?<= )\d{2}-[A-Z]{3}-\d{4}', line)
			date = checkn(date, "date", ln, '01-JAN-1979')
			
			# Division
			division = [s for s in 'BCT,PRI,ROD,MAM,VRT,INV,PLN,VRL,PHG,RNA,SYN,UNA,EST,STS,GSS,HTG,UNK'.split(',') if ' ' + s + ' ' in line]
			division = checkn(division, "division", ln, args.default_division)
			
			# Topology
			topology = [ s for s in ['linear', 'circular'] if ' ' + s + ' ' in line]
			topology = '{:<8}'.format(checkn(topology, "topology", ln, ''))
			
			# Name
			remaining = line
			for i in ['LOCUS', length, seqtype, strand, date, division, topology]:
				remaining = remaining.replace(i.strip(), '') if i is not None and not i.isspace() else remaining
			name = remaining.split()
			notname = [n for n in name if n != line.split()[1]]
			if len(notname) > 0:
				sys.exit("Line %s has unidentified value(s): %s" % (str(ln), ', '.join(notname)))
			
			name = checkn(name, "sequence name", ln, None)
			if name is None: sys.exit("Line %s has no identifiable sequence name" % (str(ln)))
			
			# Format together
				# Do name and length
			length = length.split()[0]
			llen = len(length)
			if(len(name)+llen > 27): sys.ext("Line %s sequence name %s is too long, must %s or fewer characters to allow length" % (str(ln), name, str(27-llen)))
			namelength = str('{:<' + str(28-llen) + '}').format(name) + length + ' bp '
			
				# Output
			oldline = line
			line = 'LOCUS ' + ' '*6 + namelength + strand + seqtype + ' ' + topology + ' '+ division + ' ' + date
			
			if (args.show_differences):
				sys.stdout.write("Line %s was: %s\n%s now: %s\n" % (str(ln), oldline, ' '*(5+len(str(ln))), line))
				
				for i,c in enumerate(line):
					sys.stdout.write("%s - %s\n" % (i,c))
		
		if (not args.show_differences): sys.stdout.write(line + '\n')
	
	gb.close()

exit()
