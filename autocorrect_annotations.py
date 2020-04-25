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

from Bio import SeqIO, SeqRecord, BiopythonWarning

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
	#arglist = ['-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-04-23_current/BIOD02225.gb']
	#arglist.extend("-a ATP6 -o ATP8,7 -x 50 -s N,ATG/ATA/GTG/ATT,*,C -d 53 -t 5".split(' '))
	#arglist = ['-i', '/home/thomas/MMGdatabase_currrun/1a_gbmaster_auto_run1/BIOD00100.gb', '-m', '/home/thomas/MMGdatabase_currrun/1e_nt_align/ND6.fa']
	#arglist = ['-i', '/home/thomas/MMGdatabase_testrun/1a_gbmaster_auto_run1/BIOD00100.gb', '-m', '/home/thomas/MMGdatabase_testrun/1e_nt_align/COX2.fa']
	#arglist = ['-i', '/home/thomas/MMGdatabase_testrun/1a_gbmaster_auto_run1/BIOD00100.gb', '-m', '/home/thomas/MMGdatabase_testrun/1e_nt_align/ND1.fa']
	#arglist.extend("-a ATP6 -s N,ATG/ATA/GTG/ATT,* -f N,TAG/TAA/TA,1 -d 20 -t 5 -e 1".split(' '))
	#arglist.extend("-a ATP8 -s N,ATT/ATC/AAG/ATA/TTG,* -f N,TAA/TA/T,1 -d 20 -t 5 -e 1".split(' '))
	#arglist.extend("-a COX3 -s N,ATG/ATA,* -f N,TAA/TA/T/TAG,1 -d 20 -t 5 -e 1".split(' '))
	#arglist.extend("-a ND1 -s N,TTG/ATG/ATT/TTA/ATA,* -f N,TAA/TAG/TA,1 -d 20 -t 5 -e 1".split(' '))
	#arglist.extend("-a ND4 -s N,ATG/ATA,* -f N,TAA/TA/T,1 -d 20 -t 5 -e 1".split(' '))
	#arglist.extend("-a ND6 -s N,ATT/ATA/ATC/TTG,* -f N,TAA/TAG/TA/T,1 -d 20 -t 5 -e 1".split(' '))
	#args = parser.parse_args(arglist)
	
	# Check arguments
	stringspec, overlap, alignment_distances = autocorrect_modules.check_arguments(args)
	#x = [s for s, v in alignment_distances.items() if v['body'] != v['mode']][500]
	#alignment_distances[x]
	# Read and parse gene name variants
	namevariants = autocorrect_modules.loadnamevariants()
	
	# Find universal names for inputs
	args, overlap = autocorrect_modules.standardise_names(args, overlap, namevariants)
	
	# Work through input genbank
	unrecognised_names, missing_annotation, missing_context, output_records, unidentifiable_features, context_overdist = [set(), set(), set(), [], dict(), dict()]
	
	
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
		
		if(ntf > 0):
			
			# Work through target features
			for name, feats in target_features.items():
				#name, feats = list(target_features.items())[0]
				if(args.syncronise is not None):
					autocorrect_modules.syncronise_features(name, feats, args.syncronise, seqname)
				else:
					
					feats_store = copy.deepcopy(feats)
					prev_codon_start = None
					
					for i in range(0, len(feats)):
						
						# Check if the location of this feat is identical to the prior location of the previous feat
						if(i > 0 and feats[i].location == feats_store[i-1].location):
							feats[i].location = feats[i-1].location
							if(prev_codon_start and feats[i].type == 'CDS'):
								feats[i].qualifiers['codon_start'] = prev_codon_start
							continue
						
						#i = 0
						feat = feats[i]
						
						# Set defaults
						codon_start = None
						
						# Run positional adjustments
						
						currfeat = copy.deepcopy(feat)
						posfeat = currfeat
						
						if(args.overlap is not None and ncf > 0):
							
							posfeat, record_context_overdist = autocorrect_modules.correct_positions_by_overlap(currfeat, context_features, overlap, args.maxdist, len(seq_record.seq), seqname)
							
							if(len(record_context_overdist) > 0):
								context_overdist[seqname] = record_context_overdist
							
						elif(args.match_alignment):
							
							dists = alignment_distances[seqname]
							
							# set the current feature to the correct place according to the alignment body
							posfeat, align_moved = autocorrect_modules.correct_feature_by_alignment(currfeat, stringspec, dists['body'], name, seqname, len(seq_record))
							if(dists['body'] != dists['mode']):
								posfeat = {'body': posfeat}
								align_moved = {'body': align_moved}
								posfeat['mode'], align_moved['mode'] = autocorrect_modules.correct_feature_by_alignment(currfeat, stringspec, dists['mode'], name, seqname, len(seq_record))
							
						
						currfeat = posfeat
						
						# Run string searching based adjustments
						
						if(len(stringspec) > 0):
							
							
							# Specify settings for alignment matching
								# Turn off priority for long stops 
							prioritise_long_stops = args.match_alignment is None
								# Set start reading frame based on distance if not set
							correction = [[1,3,2], [2,1,3], [3,1,2]]
							if(args.match_alignment and stringspec['start'][2] == '*' and args.force_alignment_frame):
								correction = correction[args.force_alignment_frame - 1]
								stringspec['start'][2] = correction[alignment_distances[seqname]['body']['start']%3]
							
							
							# Set up if doing multiple runs
							featdict = None
							if(type(currfeat) is dict):
								featdict = copy.deepcopy(currfeat)
								currfeat = currfeat['body']
							
							# run correct_feature_by_query
							currfeat, codon_start, query_moved = autocorrect_modules.correct_feature_by_query(currfeat, stringspec, seq_record, seqname, args.search_distance, name, args.translation_table, prioritise_long_stops)
							if(codon_start): currfeat.qualifiers['codon_start'] = codon_start
							
							# Do second run if necessary
							if(featdict is not None):
								currfeat = {'body': currfeat}
								query_moved = {'body': query_moved}
								
								
								if(stringspec['start'][2] == '*' and args.force_alignment_frame):
									correction = correction[args.force_alignment_frame - 1]
									stringspec['start'][2] = correction[alignment_distances[seqname]['mode']['start']%3]
								
								currfeat['mode'], codon_start, query_moved['mode'] = autocorrect_modules.correct_feature_by_query(featdict['mode'], stringspec, seq_record, seqname, args.search_distance, name, args.translation_table, prioritise_long_stops)
								if(codon_start): currfeat['mode'].qualifiers['codon_start'] = codon_start
								
								if(currfeat['mode'].location == currfeat['body'].location):
									currfeat = currfeat['mode']
								else:
									movement_total = {m: {e: d + query_moved[m][e] for e, d in v.items()} for m,v in align_moved.items()}
									accuracy = {m: sum([abs(alignment_distances[seqname][m][e]-d) for e, d in v.items()]) for m, v in movement_total.items()}
									bestmethod = [m for m, a in accuracy.items() if a == min(list(accuracy.values()))][0]
									currfeat = currfeat[bestmethod]
							
						
						
						# Check output and assign
						w = "Warning: new annotation for " + name + " in " + seqname
						ss = currfeat.qualifiers['codon_start']-1 if codon_start else 0
						
						if(currfeat.location.start < 0 or currfeat.location.end > len(seq_record.seq)):
							
							sys.stderr.write(w + " exceeds contig, no change made\n")
							
						elif(currfeat.location.start > currfeat.location.end):
							
							sys.stderr.write(w + " is incorrectly oriented, no change made\n")
							
						elif(autocorrect_modules.stopcount(SeqRecord.SeqRecord(currfeat.extract(seq_record.seq)[ss:]), args.translation_table, 1, False) > 0):
							
							sys.stderr.write(w + " generates internal stop codons, no change made\n")
							
						elif(currfeat.location.start != feat.location.start or 
							 currfeat.location.end != feat.location.end or 
							 ('codon_start' in feat.qualifiers and feat.qualifiers['codon_start'] != codon_start) or
							 ('codon_start' not in feat.qualifiers and codon_start is not None)):
							
							feat.location = currfeat.location
							prev_codon_start = codon_start
							
							if(codon_start and feat.type == 'CDS'):
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
