#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

import sys
import argparse
import copy
#import os
#import re

import autocorrect_modules

from Bio import SeqIO, SeqFeature, BiopythonWarning

import warnings
with warnings.catch_warnings():
	warnings.simplefilter('ignore', BiopythonWarning)

# Global variables

parser = argparse.ArgumentParser(description = "")


parser.add_argument("-i", "--input", help = "a genbank file containing one or more entries to correct", type = str, metavar = "GENBANK", required = True)

parser.add_argument("-a", "--annotation", help = "the name of the annotations to autocorrect", type = str, metavar = "ANNOT")
parser.add_argument("-o", "--overlap", help = "the context annotation and number of bases of overlap to use for correction of --annotation", type = str, metavar = "XXX,N", action = 'append')
parser.add_argument("-x", "--maxdist", help = "threshold maximum existing spacing/overlap of context/target annotations for overlap, or ungapped positions for match_alignment (default 50)", type = int, metavar = "N", default = 50)

parser.add_argument("-s", "--startstring", help = "specification of the sequence that --annotation should start with", type = str, metavar = "X,XXX[,N,X]")
parser.add_argument("-f", "--finishstring", help = "specification of the sequence that --annotation should finish with", type = str, metavar = "X,XXX[,N,X]")
parser.add_argument("-d", "--search_distance", help = "threshold maximum distance in base pairs to search for corrected start/finish", type = int, metavar = "N", default = 6)
parser.add_argument("-t", "--translation_table", help = "the amino acid translation table to use, where relevant", type = int, metavar = "N")

parser.add_argument("-y", "--syncronise", help = "the type of annotation to syncronise", type = str, metavar = "TYPE")

parser.add_argument("-m", "--match_alignment", help = "an alignment to match to when searching for a start or finish string", type = str, metavar = "ALIGN")
parser.add_argument("-e", "--force_alignment_frame", help = "force alignment matching to be guided by the supplied reading frame of the alignment", type = int, metavar = "N", choices = [1,2,3])
parser.add_argument("-w", "--show_warnings", help = "print warnings about missing or unidentifiable annotations, or annotation names that can't be recognised", action = 'store_true')

# Function definitions


# Class definitons

# Main

if __name__ == "__main__":
	
	args = parser.parse_args()
	
	# Read in arguments
	#arglist = ['-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-03-21_current/BIOD00497.gb']
	#arglist.extend("-a ND3 -o TRNG,0 -o TRNR,2 -x 35 -s N,ATA/ATT/ATG,*,C -d 45 -t 5".split(' '))
	#args = parser.parse_args(arglist)
	
	# Check arguments
	stringspec, overlap, alignment_distances, alignment_stddevs = autocorrect_modules.check_arguments(args)
	
	# Read and parse gene name variants
	namevariants = autocorrect_modules.loadnamevariants()
	
	# Find universal names for inputs
	args, overlap = autocorrect_modules.standardise_names(args, overlap, namevariants)
	
	# Work through input genbank
	unrecognised_names, missing_annotation, 	missing_context, output_records, unidentifiable_features, context_overdist = [set(), set(), set(), [], dict(), dict()]
	
	
	for seq_record in SeqIO.parse(args.input, "genbank"):
		#seq_record = next(SeqIO.parse(args.input, "genbank"))
		
		seqname = seq_record.name
		
		if(args.match_alignment is not None and seqname not in alignment_distances):
			#sys.stderr.write("Warning: no sequence for " + seqname + " can be found in " + args.match_alignment + ", this will be skipped\n")
			output_records.append(seq_record)
			continue
		
		# Get features and parse errors
		feature_names = args.annotation if args.syncronise is not None else [args.annotation]
		feature_names += overlap.keys() if overlap is not None else []
		
		features, record_unrecognised_names, record_unidentifiable_features = autocorrect_modules.get_features_from_names(seq_record, feature_names, namevariants)
		unrecognised_names.update(record_unrecognised_names)
		
		if(len(record_unidentifiable_features) > 0):
			unidentifiable_features[seqname] = record_unidentifiable_features
		
		# Assign features to categories
		target_features, context_features = [None, None]
		
		if(args.syncronise is not None):
			target_features = features
		else:
			target_features = {n:f for n, f in features.items() if n == args.annotation}
			
			if(args.overlap is not None):
				context_features = {n:f for n, f in features.items() if n in overlap.keys()}
				context_features = autocorrect_modules.check_context_features(context_features, seqname)
		
		
		# Check if have all the features needed
		ntf = len(target_features)
		ncf = sum([len(cfl) for gene, cfl in context_features.items()]) if args.overlap is not None else 1
		
		if(ntf > 0 and ncf > 0):
			
			# Work through target features
			for name, feats in target_features.items():
				#name, feats = list(target_features.items())[0]
				if(args.syncronise is not None):
					autocorrect_modules.syncronise_features(name, feats, args.syncronise, seqname)
				else:
					
					feats_store = copy.deepcopy(feats)
					
					for i in range(0, len(feats)):
						#i = 0
						
						# Check if the location of this feat is identical to the prior location of the previous feat
						if(i > 0 and feats[i].location == feats_store[i-1].location):
							feats[i].location = feats[i-1].location
							if('codon_start' in feats[i-1].qualifiers):
								feats[i].qualifiers['codon_start'] = feats[i-1].qualifiers['codon_start']
							continue
						
						feat = feats[i]
						
						# Set defaults
						codon_start = None
						
						# Run positional adjustments
						
						currfeat = copy.deepcopy(feat)
						
						if(args.overlap is not None):
							
							currfeat, record_context_overdist = autocorrect_modules.correct_positions_by_overlap(feat, context_features, overlap, args.maxdist, len(seq_record.seq), seqname)
							
							if(len(record_context_overdist) > 0):
								context_overdist[seqname] = record_context_overdist
							
						elif(args.match_alignment):
							
							# set the current feature to the correct place according to the alignment
							currfeat = autocorrect_modules.correct_feature_by_alignment(feat, stringspec, alignment_distances[seqname], name, seqname, len(seq_record))
						
						# Run string searching based adjustments
						
						if(len(stringspec) > 0):
							
							# Check whether long stops should be prioritised
							prioritise_long_stops = args.match_alignment is None
							
							# run correct_feature_by_query
							currfeat, codon_start = autocorrect_modules.correct_feature_by_query(currfeat, stringspec, seq_record, seqname, args.search_distance, name, args.translation_table, prioritise_long_stops)
						
						# Extract new start and finish
						corrected_start, corrected_finish = [currfeat.location.start, currfeat.location.end]
						
						# Check output and assign
						
						if(currfeat.location.start < 0 or currfeat.location.end > len(seq_record.seq)):
							sys.stderr.write("Warning: new annotation for " + name + " in " + seqname + " exceeds contig, no change made" + "\n")
						if(currfeat.location.start > currfeat.location.end):
							sys.stderr.write("Warning: new annotation for " + name + " in " + seqname + " is incorrectly oriented, no change made")
						elif(currfeat.location.start != feat.location.start or 
							 currfeat.location.end != feat.location.end or 
							 ('codon_start' in feat.qualifiers and feat.qualifiers['codon_start'] != codon_start)):
							feat.location = currfeat.location
							if(codon_start):
								feat.qualifiers['codon_start'] = codon_start
							
						
						# Correct truncated features
						
						if((feat.location.start == 0 or feat.location.end == len(seq_record))
							and feat.type != 'source'
							and seq_record.annotations['data_file_division'] != 'circular'):
							autocorrect_modules.correct_truncated_features(feat, len(seq_record))
							
						
					
				
			
		else:	
			if(ntf == 0): missing_annotation.add(seqname)
			if(ncf == 0): missing_context.add(seqname)
			
		
		output_records.append(seq_record)
		
	
	# Output the corrected entries
	
	if len(output_records)>0 : SeqIO.write(output_records, sys.stdout, "genbank")
	
	# Show warnings
	
	if(args.show_warnings):
		if(len(missing_annotation) > 0):
			sys.stderr.write("\nWarning, the following sequence entries were missing one or more target annotations:\n%s\n" % (', '.join(missing_annotation)))
		
		if(len(missing_context) > 0):
			sys.stderr.write("\nWarning, the following sequence entries were missing one or more context annotations:\n%s\n" % (', '.join(missing_context)))
		
		if(len(context_overdist) > 0):
			sys.stderr.write("\nWarning, the following sequence entries had context annotations that were more than " + str(args.maxdist) + " bases from target:\n")
			for seqname, cofeats in context_overdist.items():
				sys.stderr.write(seqname + ": " + ', '.join(cofeats) + "\n")
		
		if(len(unidentifiable_features) > 0):
			sys.stderr.write("\nWarning, the following sequence entries had unidentifiable annotations:\n")
			for seqname, unidfeats in unidentifiable_features.items():
				sys.stderr.write(seqname + ": " + ', '.join([f + " " + str(s) + "-" + str(e) for f, s, e in unidfeats]) + "\n")
		
		if(len(unrecognised_names) > 0):
			sys.stderr.write("\nWarning, could not recognise some feature names:\n%s\n" % (', '.join(unrecognised_names)))
	
	
	
	
	
	exit()
