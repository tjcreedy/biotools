#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

import sys
import argparse
#import os
#import re

import autocorrect_modules

from Bio import SeqIO, SeqFeature, BiopythonWarning

import warnings
with warnings.catch_warnings():
	warnings.simplefilter('ignore', BiopythonWarning)

# Global variables

parser = argparse.ArgumentParser(description = "Tool for autocorrecting the annotations in a genbank file. \nOVERLAP\nUse this option to shorten or lengthen each instance of a specified annotation according to a specified overlap with a specified 'context' annotation on a per-sequence basis. For example, ensuring that the end of annotation A is always overlapping with annotation B by 1 base. Arguments to --overlap should consist of two values, separated by a comma: the name of the annotation to overlap with, and an integer denoting the overlap desired. Positive integers imply overlap by this number of bases, a value of 0 implies exactly consecutive annotations, and negative integers imply that there are always this number of bases between the annotations. The annotation to edit is given to --annotation. Only one argument to --annotation is allowed, however multiple possible context annotations can be supplied with multiple uses of --overlap. Note that each will be applied independently so suppling multiple annotations within --overlap_maxdist that are intended to affect the same end of the --annotation may cause idiosyncratic results \nSTART/END STRING\nUse one or both of these options to shorten or lengthen each instance of a specified annotation to start or finish with a specified nucleotide or amino acid sequence. The sequence will be searched for within +/- --search_distance positions (default 10). In a comma-separated string to --startstring or --finishstring, specify the sequence type (A: Amino acid or N: Nucleotide), the sequence or list of sequences (separated by /), the codon position of the first base (N only, use * to allow any position), and a code denoting how to select which match to use if multiple matches (F, FC, C, LC, L or N). For example \"-e N,ATT,1,F\". If F is chosen, the first match will be selected starting at the position --search_distance positions before the current position and moving downstream. If FC is chosen, the first match will be selected starting at the current position and moving upstream. If C is chosen, the first match will be selected starting at the current position and moving upstream and downstream simultaneously - ties will result in no position being selected. If LC is chosen, the first match will be selected starting at the current position and moving downstream. If L is selected, the first match will be selected starting at --search_distance positions after the current position and moving upstream. Finally, if N is given, no match will be selected and the sequence output unchanged. Note: use * to represent a stop codon in AA sequences; codon position 1 is in the reading frame of the existing annotation\nSYNCRONISE\nChange the location of annotations of the type supplied to --syncronise to match the locations of any other annotations with the same name")


parser.add_argument("-i", "--input", help = "a genbank file containing one or more entries to correct", type = str, metavar = "GENBANK", required = True)

parser.add_argument("-a", "--annotation", help = "the name(s) of the annotations to autocorrect", type = str, metavar = "ANNOT", action = 'append')
parser.add_argument("-o", "--overlap", help = "the context annotation and number of bases of overlap to use for correction of --annotation", type = str, metavar = "XXX,N", action = 'append')
parser.add_argument("-x", "--overlap_maxdist", help = "threshold maximum existing spacing/overlap of context/target annotations for overlap (default 50)", type = int, metavar = "N", default = 50)

parser.add_argument("-s", "--startstring", help = "specification of the sequence that --annotation should start with", type = str, metavar = "X,XXX[,N,X]")
parser.add_argument("-f", "--finishstring", help = "specification of the sequence that --annotation should finish with", type = str, metavar = "X,XXX[,N,X]")
parser.add_argument("-d", "--search_distance", help = "threshold maximum distance in base pairs to search for corrected start/finish", type = int, metavar = "N", default = 6)
parser.add_argument("-t", "--translation_table", help = "the amino acid translation table to use, where relevant", type = int, metavar = "N")

parser.add_argument("-y", "--syncronise", help = "the type of annotation to syncronise", type = str, metavar = "TYPE")

parser.add_argument("-c", "--correct_truncated", help = "ensure annotations truncated by the end of a contig are correctly marked as partial", action = "store_true", default = False)

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
	
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-03-24_autocorrect1/BIOD00024.gb', '-a', 'COX1', '-s', 'N,ATA/ATT/ATG/ATC/ACT/ACC,*', '-m', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/nlignments_2020-03-24/COX1.fa', '-t', '5', '-e', '1', '-d', '9'])
	
	args = parser.parse_args()
	
	# Check arguments
	
	stringspec = dict()
	overlap = {}
	if(args.overlap is not None):
		if(len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --overlap")
		if(args.startstring or args.finishstring):
			sys.exit("Error: --startstring or --finishstring not compatible with --overlap. Run consecutive iterations to perform both options")
		overlap = autocorrect_modules.parse_overlap(args.overlap, args.annotation[0])
		
	elif(args.startstring is not None or args.finishstring is not None):
		if(args.annotation is None or len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --*string")
		
		stringspec = autocorrect_modules.parse_stringsearch(args.annotation[0], args.startstring, args.finishstring, args.translation_table, args.match_alignment is not None)
		
		
	elif(args.syncronise is not None):
		if(args.annotation is not None or args.overlap is not None):
			sys.exit("Error: no --annotation or --overlap should be supplied when using --syncronise")
		if(args.syncronise not in ['gene', 'CDS', 'tRNA']):
			sys.exit("Erorr: value passed to --syncronise should be gene, CDS or tRNA")
		sys.stdout.write("Running position syncronisation on %s annotations\n" % (args.syncronise))
	elif(args.correct_truncated):
		if(args.annotation is not None or args.overlap is not None):
			sys.exit("Error: no --annotation or --overlap should be supplied when using --correct_truncated")
		sys.stderr.write("Ensuring truncated annotations are properly recorded as partial\n")
	else:
		sys.exit("Error: please supply a value to --overlap, to --startstring and/or --endstring or to --syncronise, or specify --correct_truncated")
	
	# Read and parse gene name variants
	
	namevariants = autocorrect_modules.loadnamevariants()
	
	# Find universal names for inputs
	err = "Error: unrecognised locus name supplied to"
	
	if(args.syncronise is not None or args.correct_truncated):
		args.annotation = list(set(namevariants.values()))
	else:
		if(args.annotation is not None):
			if(all(a.upper() in namevariants for a in args.annotation)):
				args.annotation = [namevariants[a.upper()] for a in args.annotation]
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
	
	# Load in and process the alignment if doing matching
	if(args.match_alignment):
		alignment_distances, alignment_stddevs = autocorrect_modules.parse_alignment(args.match_alignment, args.force_alignment_frame)
	sys.stderr.write("Completed loading and parsing "+ args.match_alignment + "\n")
	
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
		
		if(args.correct_truncated and seq_record.annotations['data_file_division'] == 'circular'):
			output_records.append(seq_record)
			continue
		elif(args.match_alignment is not None and seqname not in alignment_distances):
			sys.stderr.write("Warning: no sequence for " + seqname + " can be found in " + args.match_alignment + ", this will be skipped\n")
			output_records.append(seq_record)
			continue
		
		# Get features and parse errors
		feature_names = args.annotation
		feature_names += overlap.keys() if overlap is not None else []
		
		features, record_unrecognised_names, record_unidentifiable_features = autocorrect_modules.get_features_from_names(seq_record, feature_names, namevariants)
		unrecognised_names.update(record_unrecognised_names)
		
		if(len(record_unidentifiable_features) > 0):
			unidentifiable_features[seqname] = record_unidentifiable_features
		
		# Assign features to categories
		target_features, context_features = [None, None]
		
		if(args.syncronise is not None or args.correct_truncated):
			target_features = features
		else:
			target_features = {n:f for n, f in features.items() if n == args.annotation[0]}
			
			if(args.overlap is not None):
				context_features = {n:f for n, f in features.items() if n in overlap.keys()}
				autocorrect_modules.check_context_features(context_features, args.translation_table)
		
		
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
					
					for feat in feats:
						#feat = feats[0]
						
						# Set defaults
						corrected_start, corrected_finish = [feat.location.start, feat.location.end]
						codon_start = None
						
						# Find corrections
						
						if(args.overlap is not None):
							
							corrected_start, corrected_finish, record_context_overdist = autocorrect_modules.correct_positions_by_overlap(feat, context_features, overlap, args.overlap_maxdist, len(seq_record.seq), seqname)
							
							if(len(record_context_overdist) > 0):
								context_overdist[seqname] = record_context_overdist
							
						elif(len(stringspec) > 0):
							
							if(args.match_alignment):
								
								# set the current feature to the correct place according to the alignment
								currfeat = autocorrect_modules.correct_feature_by_alignment(feat, stringspec, alignment_distances[seqname])
								# run correct_feature_by_query
								corrected_start, corrected_finish, codon_start = autocorrect_modules.correct_feature_by_query(currfeat, stringspec, seq_record, seqname, args.search_distance, name, args.translation_table, False)
							else:
								corrected_start, corrected_finish, codon_start = autocorrect_modules.correct_feature_by_query(feat, stringspec, seq_record, seqname, args.search_distance, name, args.translation_table, True)
							
						elif(args.correct_truncated):
							
							if((feat.location.start > 0 and feat.location.end < len(seq_record)) or feat.type == 'source'):
								continue
							
							corrected_start, corrected_finish = autocorrect_modules.correct_truncated_features(feat, len(seq_record))
							
						
						# Check output and assign
						if(corrected_start == feat.location.start and corrected_finish == feat.location.end):
							continue
						elif(corrected_start < 0 or corrected_finish > len(seq_record.seq)):
							sys.stderr.write("Warning: new annotation for " + name + " in " + seqname + " exceeds contig, no change made" + "\n")
							continue
						elif(corrected_start > corrected_finish):
							sys.stderr.write("Warning: new annotation for " + name + " in " + seqname + " is incorrectly oriented, no change made")
							continue
						else:
							feat.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, feat.location.strand)
							if(codon_start is not None):
								feat.qualifiers['codon_start'] = codon_start
							
						
					
				
			
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
			sys.stderr.write("\nWarning, the following sequence entries had context annotations that were more than " + str(args.overlap_maxdist) + " bases from target:\n")
			for seqname, cofeats in context_overdist.items():
				sys.stderr.write(seqname + ": " + ', '.join(cofeats) + "\n")
		
		if(len(unidentifiable_features) > 0):
			sys.stderr.write("\nWarning, the following sequence entries had unidentifiable annotations:\n")
			for seqname, unidfeats in unidentifiable_features.items():
				sys.stderr.write(seqname + ": " + ', '.join([f + " " + str(s) + "-" + str(e) for f, s, e in unidfeats]) + "\n")
		
		if(len(unrecognised_names) > 0):
			sys.stderr.write("\nWarning, could not recognise some feature names:\n%s\n" % (', '.join(unrecognised_names)))
	
	
	
	
	
	exit()
