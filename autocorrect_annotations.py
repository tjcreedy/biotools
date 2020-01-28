#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

import sys
import argparse
#import os
#import re

import urllib.request

from Bio import Seq, SeqIO, SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
from copy import copy

# Global variables

parser = argparse.ArgumentParser(description = "Tool for autocorrecting the annotations in a genbank file. \nOVERLAP\nUse this option to shorten or lengthen each instance of a specified annotation according to a specified overlap with a specified context annotation on a per-sequence basis. For example, ensuring that the end of annotation A is always overlapping with annotation B by 1 base. --overlap can be set to any integer: positive integers imply overlap by this number of bases, a value of 0 implies exactly consecutive annotations, and negative integers imply that there are always this number of bases between the annotations. The annotation to edit is given to --annotation, the context annotation (never changed) is given to --context. For --overlap, only one argument to --annotation and two arguments to --context are allowed. \nSTART/END STRING\nUse one or both of these options to shorten or lengthen each instance of a specified annotation to start or finish with a specified nucleotide or amino acid sequence. In a comma-separated string to --startstring or --finishstring, specify the sequence type (A: Amino acid or N: Nucleotide), the sequence, the reading frame of the first base (N only) and whether to use the first (F), last (L) or no (N) occurences if the sequence matches multiple times. For example \"-e N,ATT,1,F\". The sequence will be searched for within +/- --search_distance positions (default 10). Note: use * to represent a stop codon in AA sequences; reading frame 1 is the frame of the existing annotation; if the search sequence is found multiple times, the ")


parser.add_argument("-i", "--input", help = "a genbank file containing one or more entries to correct", type = str, metavar = "GENBANK", required = True)

parser.add_argument("-a", "--annotation", help = "the name(s) of the annotations to autocorrect", type = str, metavar = "ANNOT", nargs = '*')
parser.add_argument("-c", "--context", help = "the name(s) of annotations to give context to autocorrection", type = str, metavar = "CONTEXT", nargs = '*')

parser.add_argument("-o", "--overlap", help = "the number of bases of overlap to use for correction of --annotation", type = int, metavar = "N")
parser.add_argument("-m", "--overlap_maxdist", help = "threshold maximum existing spacing/overlap of context/target annotations for overlap (default 50)", type = int, metavar = "N", default = 50)

parser.add_argument("-s", "--startstring", help = "specification of the sequence that --annotation should start with", type = str, metavar = "X,XXX,N,X")
parser.add_argument("-f", "--finishstring", help = "specification of the sequence that --annotation should finish with", type = str, metavar = "X,XXX,N,X")
parser.add_argument("-d", "--search_distance", help = "threshold maximum distance in base pairs to search for corrected start/finish", type = int, metavar = "N", default = 6)

parser.add_argument("-t", "--translation_table", help = "the amino acid translation table to use, where relevant", type = int, metavar = "N")

parser.add_argument("-w", "--write_unmodified", help = "output entries that are not modified anyway", action = "store_true", default = False )
# Class definitons

# Function definitions

def loadnamevariants():
	output = {}
	for line in urllib.request.urlopen("https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt"):
		line = line.decode('utf-8').strip()
		name = line.split(";")[0]
		variants = line.split(":")[1].split(",")
		for v in variants:
			output[v] = name
	return(output)

def get_features_from_names(seqrecord, names, namevariants):
	names = [names] if isinstance(names, str) else names
	
	features = { n : list() for n in names}
	unrecognised_names = set()
	seqname = seqrecord.name
	
	for feat in seqrecord.features:
		# Remove any translations
		if('translation' in feat.qualifiers.keys()):
			del(feat.qualifiers['translation'])
		
		# Extract the feature name
		featname = 0
		if('gene' in feat.qualifiers.keys()):
			featname = feat.qualifiers['gene'][0].upper()
		elif('product' in feat.qualifiers.keys()):
			featname = feat.qualifiers['product'][0].upper()
		elif('label' in feat.qualifiers.keys()):
			featname = feat.qualifiers['label'][0].upper()
		elif(feat.type in ['source', 'misc_feature']):
			continue
		else:
			err = "Warning, %s annotation in %s" % (feat.type, seqname)
			if(hasattr(feat.location, 'start')):
				err += " from position %s to %s" % (str(int(feat.location.start)+1), str(int(feat.location.end)))
			err += " does not have a gene, product or label tag and so cannot be identified\n"
			sys.stderr.write(err)
			continue
		
		if(featname in namevariants):
			name = namevariants[featname]
			if(name in names):
				features[name].append(feat)
		else:
			unrecognised_names.add(featname)
	
	return(features, unrecognised_names)

def correct_positions(target_features, context_features, overlap, maxoverlap, seqname):
	
	# Check context_features all match in positions
	for name, feats in context_features.items():
		locations = [feat.location for feat in feats]
		if(not all(locations[0] == loc for loc in locations)):
			err = "Error, positions of " + str(len(locations)) + " context annotations for " + name + " in " + seqname + " do not match with one another:\n"
			for i, feat in enumerate(feats):
				err += "\t(" + str(i+1) + ") " + feat.type +" is located at bases " + str(int(feat.location.start)+1) + " to " + str(int(feat.location.end)) + "\n"
			sys.exit(err)
	
	# Work through combinations of target and context features
	# TODO: throw errors if annotation or context missing
	for target in target_features:
		tpos = [int(target.location.start), int(target.location.end)]
		for context_name in context_features:
			context = context_features[context_name][0]
			cpos = [int(context.location.start), int(context.location.end)]
			
			# Find all distances from the context start/end to the target start/end
			distance = list()
			for t in [0,1]:
				for c in [0,1]:
					distance.append(cpos[c]-tpos[t])
			abs_distance = [abs(d) for d in distance]
			
			# Find the meeting positions, extract the indices of the current distance and which target end this is
			min_distance_i = abs_distance.index(min(abs_distance))
			target_tpos_i = int(min_distance_i/2)
			
			# Warn and skip this context annotation if 
			if(distance[min_distance_i] > maxoverlap):
				sys.stderr.write("Warning, context annotation %s in %s is more than %s bases from target for sequence %s, this annotation will be ignored\n" % (context_name, seqname, str(maxoverlap), seqname))
				continue
			
			# Find the orientation (+ve, context follows target)
			max_distance_i = abs_distance.index(max(abs_distance))
			orientation = int(distance[max_distance_i]/abs_distance[max_distance_i])
			
			# Calculate the corrected position for the target end, from the current position plus the distance to an overlap of +1 plus the correctly-oriented overlap
			corrected_tpos = SeqFeature.ExactPosition(tpos[target_tpos_i] + distance[min_distance_i] + orientation * (overlap))
			
			# Overwrite the relevant target end position
			if(target_tpos_i == 0):
				target.location = SeqFeature.FeatureLocation(corrected_tpos, target.location.end, target.location.strand)
			else:
				target.location = SeqFeature.FeatureLocation(target.location.start, corrected_tpos, target.location.strand)

def extract_subject_region(seqrecord, feat, end, code, distance):
	'''For a given feature and end ("start" or "finish"), extract a sequence of either nucleotides (code = 'N') or amino acids (code = 'A'), consisting of the first or last position plus or minus positions equal to distance in the reading direction of the feature'''
	
	
	# Find the centre point and distances for the subject region
	central_position = feat.location.start
	distances = (distance, distance + 1)
	if((end == "start" and feat.location.strand == -1) or (end == "finish" and feat.location.strand == 1)):
		central_position = feat.location.end
		distances = (distance + 1, distance)
	
	if(code == 'A'):
		# Multiply by 3
		distances = tuple(3*d for d in distances)
		
		#Find how many trailing out-of-frame bases
		trailing_bases = len(feat.extract(seqrecord.seq)) % 3
		
		if(end == "finish"):
			modifier = 0
			if(trailing_bases is not 0):
				modifier = 3 - trailing_bases
			# Correct distances to ensure in-frame translation for out-of-frame trailing bases and for strand direction
			distances = (distances[0] + (-1 * feat.location.strand * modifier), distances[1] + (feat.location.strand * modifier))
	
	# Delimit the region and create feature
	end_positions = (central_position - distances[0], central_position + distances[1])
	subject_feat = SeqFeature(FeatureLocation(end_positions[0],end_positions[1]), strand = feat.location.strand, type = "CDS")
	
	# Extract the sequence
	sequence = subject_feat.extract(seqrecord.seq)
	
	sequence = sequence.translate(table = 5) if(code == 'A') else sequence
	
	return(sequence)

def str_is_int(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False

def find_all(a_str, sub):
	start = 0
	while True:
		start = a_str.find(sub, start)
		if start == -1: return
		yield start
		start += 1

# Testing data load
namevariants = loadnamevariants()	

seq_record_for = next(SeqIO.parse("/home/thomas/Documents/NHM_postdoc/MMGdatabase/reprocess_2019-12-18/BIOD00001.gb", "genbank"))
features_for, record_unrecognised_names = get_features_from_names(seq_record_for, "ATP8", namevariants)
features_for = features_for["ATP8"]

seq_record_rev = next(SeqIO.parse("/home/thomas/Documents/NHM_postdoc/MMGdatabase/reprocess_2019-12-18/BIOD00002.gb", "genbank"))
features_rev, record_unrecognised_names = get_features_from_names(seq_record_rev, "ATP8", namevariants)
features_rev = features_rev["ATP8"]

feat = copy(features_for[0])
feat = copy(features_rev[0])

feat.location = FeatureLocation(feat.location.start, feat.location.end+1, strand = 1)
feat.location = FeatureLocation(feat.location.start, feat.location.end+2, strand = 1)

feat.location = FeatureLocation(feat.location.start-1, feat.location.end, strand = -1)
feat.location = FeatureLocation(feat.location.start-2, feat.location.end, strand = -1)

distance = 3
seqrecord = seq_record_for
seqrecord = seq_record_rev
end = "start"
end = "finish"
code = 'N'
code = 'A'
				



if __name__ == "__main__":
	
	# Read in arguments
	#args = parser.parse_args(['-a', "NAD2", '-c', 'TRNM(CAU)', '-o', 0, '-i', 'source/BIOD00005.gb', '-w'])
	args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/reprocess_2019-12-18/BIOD00001.gb', '-a', 'ND5', '-c' 'TRNM(CAU)', '-t', '5'])
	args = parser.parse_args()
	
	# Check arguments
	
	stringspec = dict()
	
	if(args.overlap is not None):
		if(len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --overlap")
		if(len(args.context) not in [1,2]):
			sys.exit("Error: please supply one or two context annotation names to --context if using --overlap")
		if(args.startstring or args.finishstring):
			sys.exit("Error: --startstring or --finishstring not compatible with --overlap. Run consecutive iterations to perform both options")
	if(args.startstring is not None or args.finishstring is not None):
		for end, searchstring in zip(['start','finish'], [args.startstring, args.finishstring]):
			if(searchstring is None):
				continue
			search = tuple(searchstring.split(","))
			err = "Error: search string " + searchstring
			if(search[0] not in ["A", "N"]):
				sys.exit(err + " has an unrecognised sequence type")
			elif(search[0] == 'A' and len(search) not in [2,3]):
				sys.exit(err + " does not have only two or three items")
			elif(search[0] == 'N' and len(search) not in [3,4]):
				sys.exit(err + " does not have only three or four items")
			elif(search[0] == "A" and search[1] not in list("GPAVLIMCFYWHKRQNEDST*")):
				sys.exit(err + " is specified as amino acid but non-standard character(s) included (must be GPAVLIMCFYWHKRQNEDST*)")
			elif(search[0] == "N" and search[1] not in list("ATGC")):
				sys.exit(err + " is specified as nucleotide but non-standard character(s) included (must be ATGC)")
			elif(search[0] == "N" and (not str_is_int(search[2]) or search[2] not in ['1','2','3'])):
				sys.exit(err + " has an unrecognised reading frame")
			elif(search[0] == "A" and str_is_int(search[2])):
				sys.exit(err + " is searching for an amino acid sequence but includes what seems to be a reading frame")
			elif(( (search[0] == "N" and len(search) == 4) or (search[0] == "A" and len(search) == 3) ) and (not search[3] in ['F','L','N'])):
				sys.exit(err + " has an unrecognised multiple-hit handling instruction ( must be F, L or N)")
			stringspec[end] = search
	else:
		sys.exit("Error: please supply a value to --overlap or to --startstring and/or --endstring")
	
	# Read and parse gene name variants
	
	namevariants = loadnamevariants()	
	
	# Find universal names for inputs
	err = "Error: unrecognised locus name supplied to"
	if all(a.upper() in namevariants for a in args.annotation):
		args.annotation = [namevariants[a.upper()] for a in args.annotation]
	else:
		err = err + " --annotation"
		sys.exit(err)
	
	if(args.overlap is not None):
		if all(c.upper() in namevariants for c in args.context):
			args.context = [namevariants[c.upper()] for c in args.context]
		else:
			#die with error
			err = err + " --context"
			sys.exit(err)
	
	# Work through input genbank
	
	unrecognised_names = set()
	output_records = list()
	
	for seq_record in SeqIO.parse(args.input, "genbank"):
		#seq_record = next(SeqIO.parse(args.input, "genbank"))
		seqname = seq_record.name
		write = False
		features, record_unrecognised_names = get_features_from_names(seq_record, args.annotation + args.context, namevariants)
		unrecognised_names.update(record_unrecognised_names)
		
		if(args.overlap is not None):
			target_features = features.pop(args.annotation[0])
			context_features = features
			
			ntf = len(target_features)
			ncf = sum([len(cfl) for gene, cfl in context_features.items()])
			if(ntf > 0 and ncf > 0):
				correct_positions(target_features, context_features, args.overlap, args.overlap_maxdist, seqname)
				write = True
			else:
				err = "Warning, sequence " + seqname + " has"
				if ntf == 0: err += " no " + args.annotation[0] + " annotation" 
				if ntf == 0 and ncf == 0: err += " and"
				if ncf == 0: err += " none of the specified context annotation(s)"
				sys.stderr.write(err + "\n")
				write = args.write_unmodified
		
		elif(len(stringspec) > 0):
			features = features[args.annotation[0]]
			
			nf = len(features)
			if(nf > 0):
				feat = features[0]
				
				# Store original start and end
				instart, infin = feat.location.start, feat.location.end
				
				# Search for changes
				startchange, finchange = 0, 0
				
				errstart = "Warning: sequence " + seqname
				
				for end, search in stringspec.items():
					end = "finish"
					search = ('A','*','F')
					
					# Check if already ends with the searched sequence
					nuc_seq = feat.extract(seq_record.seq)
					if(search[0] == 'N'):
						if((end is "start" and search[2] == 1 and nuc_seq.startswith(search[1])) or  # Starts with sequence, rf is 1
						   (end is "finish" and nuc_seq.endswith(search[1]) and nuc_seq.rfind(search[1]) % 3 + 1 == search[2])): # ends with and is in the appropriate frame
							continue
					else:
						aa_seq = nuc_seq.translate(table = args.translation_table)
						if((end is "start" and aa_seq.startswith(search[1])) or
						    end is "finish" and aa_seq.endswith(search[1])):
							continue
					
					# Generate sequence for searching
					subject_sequence = extract_subject_region(seq_record, feat, end, search[0], args.search_distance)
					
					
					# TODO FROM HERE
					
					# Find locations of query
					locations = list(find_all(str(subject_sequence), search[1]))
					
					# Check if in reading frame
					if(search[0] == 'N'):
						locations = [l for l in locations if (l + 1) % 3 == int(search[2]) + (args.search_distance % 3)]
						# Retain location if that location's rf (l+1)%3 is equal to the (target rf converted to subject rf)
					
					# Parse location results
					errend = search[1]  + " at the " + end + " of " + args.annotation[0] + "\n"
					if(len(locations) > 0):
						location = locations[0]
						if(len(locations) > 1):
							if(start[3] == 'N'):
								sys.stderr.write(errstart + " has multiple matches of " + errend)
							elif(start[3] == 'L'):
								location = locations[-1]
						if(search[0] == 'A'):
							
							if(end is "start):
								startchange = location - changemod / changemult
							else:
								finchange = location - changemod / changemult
					else:
						sys.stderr.write(errstart + " has no matches of " + errend)
					
				
				if(startchange == 0 and finchange == 0):
					write = args.write_unmodified
				
				feat.location = SeqFeature.FeatureLocation(
						SeqFeature.ExactPosition(instart + startchange),
						SeqFeature.ExactPosition(infin + startchange),
						feat.location.strand)
				
				
				
				
				
				
				
				
				
				
				.translate(table=5)
				
				
			else:
				err = "Warning, sequence " + seqname + " has no " + args.annotation[0] + " annotation(s)\n"
				sys.stderr.write(err)
				write = args.write_unmodified
		
		
		
		
		if write: output_records.append(seq_record) 
	
	if len(output_records)>0 : SeqIO.write(output_records, sys.stdout, "genbank")
	
	if(len(unrecognised_names) > 0):
		sys.stderr.write("Warning, could not recognise some feature names - %s \n" % (', '.join(unrecognised_names)))
	
	exit()