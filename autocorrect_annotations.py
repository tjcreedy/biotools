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
	
	# Read in arguments
	
	#args = parser.parse_args(['-a', "ATP6", '-c', 'ATP8', '-o', '7', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-04_2edited/SPSO00168.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "CYTB", '-c', 'TRNS(UGA)', '-o', '2', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00109.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "NAD2", '-o', 'TRNW,2', '-o', 'TRNS,-20', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00001.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "NAD", '-c', 'TRNH(GUG)', '-o', '0', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00409.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "NAD5", '-o', 'TRNF,1', '-m', '100', '-t', '5', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00471.gb'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/QINL005.gb', '-a', 'ND5', '-s', 'N,ATT/ATA/ATG/ATC,1,L', '-d', '3', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/testing/BIOD00550.gb', '-a', 'ND2', '-s', 'N,ATA/ATG/ATC/TTG/ATT,*,LC', '-d', '20', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/CCCP00094.gb', '-a', 'NAD1', '-f', 'N,TAA/TAG,1,F', '-d', '220', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-04_2edited/BIOD00622.gb', '-a', 'ATP6', '-f', 'N,TAG/TAA/TA,1,F', '-d', '15', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00010.gb', '-a', 'NAD4', '-s', 'N,ATG/ATA,1,F', '-d', '6', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/CCCP00017.gb', '-y', 'gene'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00109.gb', '-a', 'ND6', '-s', 'N,ATT/ATA/ATC/TTG/TTT,*,FC', '-d', '6', '-t', '5'])
	
	#args = parser.parse_args(['-a', "NAD2", '-o', 'TRNM,3', '-o', 'TRNI,-18', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00002.gb', '-m', '50', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/testin/BIOD01170.gb', '-a', 'NAD2', '-s', 'N,ATA/ATG/ATC/TTG/ATT,*,LC', '-d', '20', '-t', '5'])
	#args = parser.parse_args(['-a', "NAD2", '-o', 'TRNW,2', '-o', 'TRNS,2', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/testin/BIOD01170.gb', '-m', '50', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/testin/BIOD01170.gb', '-a', 'NAD2', '-f', 'N,TAA/TA,1,F, '-d', '21', '-t', '5'])
	
	#args = parser.parse_args(['-a', "ATP6", '-o', 'ATP8,7', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/testin/BIOD01796.gb', '-m', '50'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/testin/BIOD01796.gb', '-a', 'ATP6', '-s', 'N,TTG/ATG/ATA,*,C', '-d', '53', '-t', '5'])
	
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-03-21_current/GBDL01179.gb', '-a', 'COX3', '-o', 'ATP6,1', '-x', '50', '-s', 'N,ATG/ATA,*,C', '-d', '28', '-t', '5'])
	
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/testing/gbmaster/BIOD00821.gb', '-d', '9', '-e', '1', '-m', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/testing/nt_align_reduced/ND5.fa', '-a', 'ND5', '-s', 'N,ATT/ATA/ATG/ATC,1', '-f', 'N,TAA/TA/T,1', '-d', '30', '-t', '5'])
	
	args = parser.parse_args()
	
	# Check arguments
	
	stringspec = dict()
	overlap = {}
	
	if(args.syncronise is not None):
		if(args.startstring is not None or
			args.finishstring is not None or
			args.overlap is not None or 
			args.match_alignment is not None):
			sys.exit("Error: --syncronise is not compatible with --startstring, --finishstring, --annotation, --overlap or --match_alignment")
		
		if(args.syncronise not in ['gene', 'CDS', 'tRNA']):
			sys.exit("Error: value passed to --syncronise should be gene, CDS or tRNA")
			sys.stdout.write("Running position syncronisation on %s annotations\n" % (args.syncronise))
		
	elif(args.annotation is not None):
		
		if(args.overlap is not None):
			
			if(args.match_alignment is None):
				overlap = autocorrect_modules.parse_overlap(args.overlap, args.annotation)
			else:
				sys.exit("Error: --match_alignment and --overlap are mutually exclusive")
			
		elif(args.match_alignment is not None):
			
			# Load in and process the alignment if doing matching
			alignment_distances, alignment_stddevs = autocorrect_modules.parse_alignment(args.match_alignment, args.force_alignment_frame)
			sys.stderr.write("Completed loading and parsing "+ args.match_alignment + "\n")
			
		elif(args.startstring is None and args.finishstring is None):
			
			sys.exti("Error: insufficient arguments - are you missing at least one of --overlap, --match_alignment, --startstring and/or --finishstring?")
		
		
		if(args.startstring is not None or args.finishstring is not None):
			stringspec = autocorrect_modules.parse_stringsearch(args.annotation, args.startstring, args.finishstring, args.translation_table, args.match_alignment is not None)
		
	else:
		sys.exit("Error: insufficient arguments - are you missing --annotation?")
	
	
	# Read and parse gene name variants
	
	namevariants = autocorrect_modules.loadnamevariants()
	
	# Find universal names for inputs
	err = "Error: unrecognised locus name supplied to"
	
	if(args.syncronise is not None):
		args.annotation = list(set(namevariants.values()))
	else:
		if(args.annotation is not None):
			if(args.annotation.upper() in namevariants):
				args.annotation = namevariants[args.annotation.upper()]
			else:
				err = err + " --annotation"
				sys.exit(err)
		
		if(overlap is not None):
			if all(c.upper() in namevariants for c in overlap.keys()):
				overlap = { namevariants[c.upper()]:o for c,o in overlap.items() }
			else:
				#die with error
				err = err + " --overlap"
				sys.exit(err)
	
	
	# Work through input genbank
	
	unrecognised_names = set()
	missing_annotation = set()
	missing_context = set()
	output_records = list()
	unidentifiable_features = dict()
	context_overdist = dict()
	
	
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
				autocorrect_modules.check_context_features(context_features, seqname)
		
		
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
							continue
						
						feat = feats[i]
						
						# Set defaults
						corrected_start, corrected_finish = [feat.location.start, feat.location.end]
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
							prioritise_long_stops = args.match_alignment is not None
							
							# run correct_feature_by_query
							currfeat, codon_start = autocorrect_modules.correct_feature_by_query(currfeat, stringspec, seq_record, seqname, args.search_distance, name, args.translation_table, prioritise_long_stops)
						
						# Extract new start and finish
						corrected_start, corrected_finish = [currfeat.location.start, currfeat.location.end]
						
						# Check output and assign
						
						if(corrected_start < 0 or corrected_finish > len(seq_record.seq)):
							sys.stderr.write("Warning: new annotation for " + name + " in " + seqname + " exceeds contig, no change made" + "\n")
						if(corrected_start > corrected_finish):
							sys.stderr.write("Warning: new annotation for " + name + " in " + seqname + " is incorrectly oriented, no change made")
						elif(corrected_start != feat.location.start or corrected_finish != feat.location.end):
							feat.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, feat.location.strand)
							if(codon_start is not None):
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
