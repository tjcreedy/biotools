#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

import sys
import argparse
#import os
#import re

import urllib.request

from Bio import SeqIO, SeqFeature
from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
	warnings.simplefilter('ignore', BiopythonWarning)

# Global variables

parser = argparse.ArgumentParser(description = "Tool for autocorrecting the annotations in a genbank file. \nOVERLAP\nUse this option to shorten or lengthen each instance of a specified annotation according to a specified overlap with a specified context annotation on a per-sequence basis. For example, ensuring that the end of annotation A is always overlapping with annotation B by 1 base. --overlap can be set to any integer: positive integers imply overlap by this number of bases, a value of 0 implies exactly consecutive annotations, and negative integers imply that there are always this number of bases between the annotations. The annotation to edit is given to --annotation, the context annotation (never changed) is given to --context. For --overlap, only one argument to --annotation and two arguments to --context are allowed. \nSTART/END STRING\nUse one or both of these options to shorten or lengthen each instance of a specified annotation to start or finish with a specified nucleotide or amino acid sequence. The sequence will be searched for within +/- --search_distance positions (default 10). In a comma-separated string to --startstring or --finishstring, specify the sequence type (A: Amino acid or N: Nucleotide), the sequence, the codon position of the first base (N only), and a code denoting how to select which match to use if multiple matches (F, FC, C, LC, L or N). For example \"-e N,ATT,1,F\". If F is chosen, the first match will be selected starting at the position --search_distance positions before the current position and moving downstream. If FC is chosen, the first match will be selected starting at the current position and moving upstream. If C is chosen, the first match will be selected starting at the current position and moving upstream and downstream simultaneously - ties will result in no position being selected. If LC is chosen, the first match will be selected starting at the current position and moving downstream. If L is selected, the first match will be selected starting at --search_distance positions after the current position and moving upstream. Finally, if N is given, no match will be selected and the sequence output unchanged. Note: use * to represent a stop codon in AA sequences; codon position 1 is in the reading frame of the existing annotation")


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

def correct_positions_by_overlap(target_features, context_features, overlap, maxoverlap, seqname):
	
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

def correct_feature_by_query(feat, query_spec, seq_record, seqname, distance, featurename):
	#query_spec = stringspec
	#distance = args.search_distance
	#featurename = args.annotation[0]
	
	feat_start, feat_finish = feat.location.start, feat.location.end
	errstart = "Warning: sequence " + seqname
	
	for end, search in query_spec.items():
		#end, search = next(iter(query_spec.items()))
		
		# Unpack search tuple
		code, query, out_rf, selector = search if len(search) == 4 else search + tuple("X")
		selector = out_rf if code == 'A' else selector
		
		# Check if already ends with the searched sequence
		if(end_already_correct(feat.extract(seq_record.seq), query, end, code, out_rf)):
			continue
		
		# Generate sequence for searching
		subject_sequence, subject_start = extract_subject_region(seq_record, feat, end, code, distance)
		
		# Find locations of query
		locations = list(find_all(str(subject_sequence), query))
		
		# Retain only locations in the specified reading frame
		if(code == 'N'):
			if(end == "start"):
				locations = [l for l in locations if (l + distance) % 3 + 1 == int(out_rf)]
				# Retain location if that location's rf (l+1)%3 is equal to the (target rf converted to subject rf)
			else:
				locations = [l for l in locations if (len(feat)-distance + l - 1) % 3 + 1 == int(search[2])]
		
		# Parse location results
		errend = query  + " at the " + end + " of " + featurename + "\n"
		
		if(len(locations) > 0):
			# Select the first location by default
			location = locations[0]
			
			if(len(locations) > 1):
				if(selector == 'N'):
					# If more than 1 locations but the user wishes none to be selected
					sys.stderr.write(errstart + " has multiple matches of " + errend)
					location = distance
				elif(selector == 'L'):
					# If more than 1 locations but the user wants the last selected
					location = locations[-1]
				elif(selector == 'C'):
					loc_dist = [abs(distance-l) for l in locations]
					if(loc_dist.count(min(loc_dist)) == 1):
						location = locations[loc_dist.index(min(loc_dist))]
					else:
						sys.stderr.write(errstart + "has multiple closest matches of " + errend)
						location = distance
				elif(selector == "FC"):
					locations = [l for l in locations if l < distance]
					if(len(locations) > 0):
						location = locations[-1]
					else:
						sys.stderr.write(errstart + "has no first closest matches of " + errend)
						location = distance
				elif(selector == "LC"):
					locations = [l for l in locations if l > distance]
					if(len(locations) > 0):
						location = locations[0]
					else:
						sys.stderr.write(errstart + "has no last closest matches of " + errend)
						location = distance
			
			# Convert locations if on reverse strand
			if(feat.location.strand == -1):
				locations = [abs(l - (len(subject_sequence))) for l in locations]
			
			# Generate the new end position
				# Correct by one base if at the finish end
			change = location + feat.location.strand if end == "finish" else location
				# Multiply by 3 if AA
			change = change * 3 if code == 'A' else change
				# Calculate
			newend = subject_start + change
			
			if((end == "start" and feat.location.strand == 1) or (
					end == "finish" and feat.location.strand == -1)):
				feat_start = newend
			else:
				feat_finish = newend
			
		else:
			sys.stderr.write(errstart + " has no matches of " + errend)
	
	return(feat_start, feat_finish)

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
			if(trailing_bases != 0):
				modifier = 3 - trailing_bases
			# Correct distances to ensure in-frame translation for out-of-frame trailing bases and for strand direction
			distances = (distances[0] + (-1 * feat.location.strand * modifier), distances[1] + (feat.location.strand * modifier))
	
	# Delimit the region and create feature
	end_positions = (central_position - distances[0], central_position + distances[1])
	subject_feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(end_positions[0],end_positions[1]), strand = feat.location.strand, type = "CDS")
	
	# Extract the sequence
	sequence = subject_feat.extract(seqrecord.seq)
	
	sequence = sequence.translate(table = 5) if(code == 'A') else sequence
	
	return(sequence, end_positions[0])

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

def end_already_correct(nuc_seq, query_seq, end, code, frame):
	#nuc_seq = feat.extract(seq_record.seq)
	#query_seq = query
	#frame = out_rf
	
	
	if(code == 'N'):
		return((end == "start" and frame == "1" and nuc_seq.startswith(query_seq)) or # Starts with sequence, rf is 1
			(end == "finish" and nuc_seq.endswith(query_seq) and nuc_seq.rfind(query_seq) % 3 + 1 == int(frame)))
	else:
		aa_seq = nuc_seq.translate(table = args.translation_table)
		return((end == "start" and aa_seq.startswith(query_seq)) or
			(end == "finish" and aa_seq.endswith(query_seq)))

if __name__ == "__main__":
	
	# Read in arguments
	#args = parser.parse_args(['-a', "NAD2", '-c', 'TRNM(CAU)', '-o', 0, '-i', 'source/BIOD00005.gb', '-w'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/reprocess_2019-12-18/BIOD00001.gb', '-a', 'ND5', '-s', 'A,M,F', '-f', 'N,TAG,1,C', '-t', '5'])
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
		if(args.annotation is None or len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --*string")
		if(args.context is not None):
			sys.exit("Error: no context annotations are used in --*string")
		
		for end, searchstring in zip(['start','finish'], [args.startstring, args.finishstring]):
			if(searchstring is None):
				continue
			search = tuple(searchstring.split(","))
			err = "Error: search string " + searchstring
			if(search[0] not in ['A', 'N']):
				sys.exit(err + " has an unrecognised sequence type")
			elif(search[0] == 'A'):
				if(len(search) != 3):
					sys.exit(err + " does not have only two or three items")
				elif(any(s not in list("GPAVLIMCFYWHKRQNEDST*") for s in search[1])):
					sys.exit(err + " is specified as amino acid but non-standard character(s) included (must be GPAVLIMCFYWHKRQNEDST*)")
				elif(str_is_int(search[2])):
					sys.exit(err + " is searching for an amino acid sequence but includes what seems to be a reading frame")
				elif(args.translation_table is None):
					sys.exit(err + " is specified as amino acid but no --translation_table given")
			else:
				if(len(search) != 4):
					sys.exit(err + " does not have only three or four items")
				elif(any(s not in list("ATGC") for s in search[1])):
					sys.exit(err + " is specified as nucleotide but non-standard character(s) included (must be ATGC)")
				elif(not str_is_int(search[2]) or search[2] not in ['1','2','3']):
					sys.exit(err + " has an unrecognised reading frame")
			
			mmhi = search[3] if search[0] == 'N' else search[2]
			if( mmhi not in ['F','L','C','FC','LC','N'] ):
				sys.exit(err + " has an unrecognised multiple-hit handling instruction ( must be F, FC, C, LC, L or N)")
			
			stringspec[end] = search
		
		args.context = []
		
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
				correct_positions_by_overlap(target_features, context_features, args.overlap, args.overlap_maxdist, seqname)
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
				for feat in features:
					#feat = features[0]
					
					# Get new start and finish positions
					corrected_start, corrected_finish = correct_feature_by_query(feat, stringspec, seq_record, seqname, args.search_distance, args.annotation[0])
					
					if(corrected_start == feat.location.start and corrected_finish == feat.location.end):
						write = args.write_unmodified
						continue
					else:
						feat.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, feat.location.strand)
						write = True
			else:
				err = "Warning, sequence " + seqname + " has no " + args.annotation[0] + " annotation(s)\n"
				sys.stderr.write(err)
				write = args.write_unmodified
		
		if write: output_records.append(seq_record) 
	
	if len(output_records)>0 : SeqIO.write(output_records, sys.stdout, "genbank")
	
	if(len(unrecognised_names) > 0):
		sys.stderr.write("Warning, could not recognise some feature names - %s \n" % (', '.join(unrecognised_names)))
	
	exit()