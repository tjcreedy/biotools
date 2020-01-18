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

# Global variables

parser = argparse.ArgumentParser(description = "Tool for autocorrecting the annotations in a genbank file. \nOVERLAP\nCurrently, the only option available is to shorten or lengthen each instance of a specified annotation according to a specified overlap with a specified context annotation on a per-sequence basis. For example, ensuring that the end of annotation A is always overlapping with annotation B by 1 base. --overlap can be set to any integer: positive integers imply overlap by this number of bases, a value of 0 implies exactly consecutive annotations, and negative integers imply that there are always this number of bases between the annotations. The annotation to edit is given to --annotation, the context annotation (never changed) is given to --context. For --overlap, only one argument to --annotation and two arguments to --context are allowed.")


parser.add_argument("-i", "--input", help = "a genbank file containing one or more entries to correct", type = str, metavar = "GENBANK", required = True)

parser.add_argument("-a", "--annotation", help = "the name(s) of the annotations to autocorrect", type = str, metavar = "ANNOT", nargs = '*')
parser.add_argument("-c", "--context", help = "the name(s) of annotations to give context to autocorrection", type = str, metavar = "CONTEXT", nargs = '*')

parser.add_argument("-o", "--overlap", help = "the number of bases of overlap to use for correction of --annotation", type = int, metavar = "N")
parser.add_argument("-m", "--overlap_maxdist", help = "threshold maximum existing spacing/overlap of context/target annotations for overlap (default 50)", type = int, metavar = "N", default = 50)

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

def get_features_from_names(seqrecord, target_names, context_names, namevariants):
	target_features = list()
	context_features = { c : list() for c in context_names}
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
			sys.stderr.write("Warning, %s annotation in %s from position %s to %s does not have a gene, product or label tag and so cannot be identified\n" % (feat.type, seqname, str(int(feat.location.start)+1), str(int(feat.location.end))))
			continue
		
		if(featname in namevariants):
			name = namevariants[featname]
			if(name in target_names):
				target_features.append(feat)
			elif(name in context_names):
				context_features[name].append(feat)
		else:
			unrecognised_names.add(featname)
	
	return(target_features, context_features, unrecognised_names)

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

if __name__ == "__main__":
	
	# Read in arguments
	#args = parser.parse_args(['-a', "NAD2", '-c', 'TRNM(CAU)', '-o', 0, '-i', 'source/BIOD00005.gb', '-w'])
	args = parser.parse_args()
	
	# Check arguments
	
	if(args.overlap != None):
		if(len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --overlap")
		if(len(args.context) not in [1,2]):
			sys.exit("Error: please supply one or two context annotation names to --contex if using --overlap")
	else:
		sys.exit("Error: please supply a value to --overlap")
	
	# Read and parse gene name variants
	
	namevariants = loadnamevariants()	
	
	# Find universal names for inputs
	err = "Error: unrecognised locus name supplied to"
	if all(a.upper() in namevariants for a in args.annotation):
		args.annotation = [namevariants[a.upper()] for a in args.annotation]
	else:
		err = err + " --annotation"
		sys.exit(err)
	
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
		target_features, context_features, record_unrecognised_names = get_features_from_names(seq_record, args.annotation, args.context, namevariants)
		unrecognised_names.update(record_unrecognised_names)
		
		if(args.overlap != None):
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
		
		if write: output_records.append(seq_record) 
	
	if len(output_records)>0 : SeqIO.write(output_records, sys.stdout, "genbank")
	
	if(len(unrecognised_names) > 0):
		sys.stderr.write("Warning, could not recognise some feature names - %s \n" % (', '.join(unrecognised_names)))
	
	exit()