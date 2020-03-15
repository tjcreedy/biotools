#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

# Imports

import sys
import argparse
#import os
#import re
import urllib.request


from collections import defaultdict
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
	warnings.simplefilter('ignore', BiopythonWarning)

# Global variables

parser = argparse.ArgumentParser(description = "Tool for autocorrecting the annotations in a genbank file. \nOVERLAP\nUse this option to shorten or lengthen each instance of a specified annotation according to a specified overlap with a specified 'context' annotation on a per-sequence basis. For example, ensuring that the end of annotation A is always overlapping with annotation B by 1 base. Arguments to --overlap should consist of two values, separated by a comma: the name of the annotation to overlap with, and an integer denoting the overlap desired. Positive integers imply overlap by this number of bases, a value of 0 implies exactly consecutive annotations, and negative integers imply that there are always this number of bases between the annotations. The annotation to edit is given to --annotation. Only one argument to --annotation is allowed, however multiple possible context annotations can be supplied with multiple uses of --overlap. Note that each will be applied independently so suppling multiple annotations within --overlap_maxdist that are intended to affect the same end of the --annotation may cause idiosyncratic results \nSTART/END STRING\nUse one or both of these options to shorten or lengthen each instance of a specified annotation to start or finish with a specified nucleotide or amino acid sequence. The sequence will be searched for within +/- --search_distance positions (default 10). In a comma-separated string to --startstring or --finishstring, specify the sequence type (A: Amino acid or N: Nucleotide), the sequence or list of sequences (separated by /), the codon position of the first base (N only, use * to allow any position), and a code denoting how to select which match to use if multiple matches (F, FC, C, LC, L or N). For example \"-e N,ATT,1,F\". If F is chosen, the first match will be selected starting at the position --search_distance positions before the current position and moving downstream. If FC is chosen, the first match will be selected starting at the current position and moving upstream. If C is chosen, the first match will be selected starting at the current position and moving upstream and downstream simultaneously - ties will result in no position being selected. If LC is chosen, the first match will be selected starting at the current position and moving downstream. If L is selected, the first match will be selected starting at --search_distance positions after the current position and moving upstream. Finally, if N is given, no match will be selected and the sequence output unchanged. Note: use * to represent a stop codon in AA sequences; codon position 1 is in the reading frame of the existing annotation\nSYNCRONISE\nChange the location of annotations of the type supplied to --syncronise to match the locations of any other annotations with the same name")


parser.add_argument("-i", "--input", help = "a genbank file containing one or more entries to correct", type = str, metavar = "GENBANK", required = True)

parser.add_argument("-a", "--annotation", help = "the name(s) of the annotations to autocorrect", type = str, metavar = "ANNOT", action = 'append')
parser.add_argument("-o", "--overlap", help = "the context annotation and number of bases of overlap to use for correction of --annotation", type = str, metavar = "XXX,N", action = 'append')
parser.add_argument("-m", "--overlap_maxdist", help = "threshold maximum existing spacing/overlap of context/target annotations for overlap (default 50)", type = int, metavar = "N", default = 50)

parser.add_argument("-s", "--startstring", help = "specification of the sequence that --annotation should start with", type = str, metavar = "X,XXX,N,X")
parser.add_argument("-f", "--finishstring", help = "specification of the sequence that --annotation should finish with", type = str, metavar = "X,XXX,N,X")
parser.add_argument("-d", "--search_distance", help = "threshold maximum distance in base pairs to search for corrected start/finish", type = int, metavar = "N", default = 6)
parser.add_argument("-t", "--translation_table", help = "the amino acid translation table to use, where relevant", type = int, metavar = "N")

parser.add_argument("-y", "--syncronise", help = "the type of annotation to syncronise", type = str, metavar = "TYPE")

parser.add_argument("-w", "--show_warnings", help = "print warnings about missing or unidentifiable annotations, or annotation names that can't be recognised", action = 'store_true')

# Class definitons

# Function definitions

def loadnamevariants():
	output = {}
	for line in urllib.request.urlopen("https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt"):
		line = line.decode('utf-8').strip()
		name = line.split(";")[0]
		annotype = line.split(":")[0].split(";")[1]
		variants = line.split(":")[1].split(",")
		for v in variants:
			for g in ['', ' ']:
				v = v.replace(g, '')
				for s in ['',' GENE', ' '+annotype.upper()]:
					output[v+s] = name
	return(output)

def get_features_from_names(seqrecord, names, namevariants):
	#names = args.annotation + args.context or []
	#seqrecord = seq_record
	names = [names] if isinstance(names, str) else names
	
	features = defaultdict(list)
	unrecognised_names = set()
	unidentifiable_features = set()
	#seqname = seqrecord.name
	
	for feat in seqrecord.features:
		#feat = seqrecord.features[0]
		# Remove any translations
		if('translation' in feat.qualifiers.keys()):
			del(feat.qualifiers['translation'])
		
		# Extract the feature name
		featname = 0
		nametags = ['gene', 'product', 'label', 'standard_name']
		if(any(t in feat.qualifiers.keys() for t in nametags)):
			for t in nametags:
				if(t in feat.qualifiers.keys()):
					featname = feat.qualifiers[t][0].upper()
					break
		elif(feat.type in ['source', 'misc_feature', 'repeat_region', 'D-loop', 'rep_origin','gap']):
			continue
		else:
			unidentifiable_features.add((feat.type, feat.location.start, feat.location.end))
			#err = "Warning, can't identify %s %s annotation" % (seqname, feat.type)
			#if(hasattr(feat.location, 'start')):
			#	err += " %s-%s" % (str(int(feat.location.start)+1), str(int(feat.location.end)))
			#err += ": no gene/product/label/standard_name tag\n"
			#sys.stderr.write(err)
			continue
		
		if(featname in namevariants):
			name = namevariants[featname]
			if(name in names):
				features[name].append(feat)
		else:
			unrecognised_names.add(featname)
	
	return(features, unrecognised_names, unidentifiable_features)

def syncronise_features(features, synctype, seqname):
	#synctype = 'gene'
	for name, feats in features.items():
		#name, feats = list(features.items())[2]
		# Organise features
		target_feats = []
		other_feats = []
		
		for feat in feats:
			if(feat.type == synctype):
				target_feats.append(feat)
			else:
				other_feats.append(feat)
		
		warnstart = "Warning, " + seqname + " does not have any "
		warnend = " annotations of " + name + "\n"
		if(len(target_feats) < 1):
			#sys.stderr.write(warnstart + synctype + warnend)
			continue
		if(len(other_feats) < 1):
			#sys.stderr.write(warnstart + "non-" + synctype + warnend)
			continue
		
		# Check that the other features are all the same
		if(not all(other_feats[0].location == feat.location for feat in other_feats)):
			sys.stderr.write("Warning, positions of " + str(len(other_feats)) + " non-" + synctype + " annotations for " + name + " in " + seqname + " do not match, this entry will not be modified\n")
			continue
		
		# Correct the target features
		for feat in target_feats:
			feat.location = other_feats[0].location

def correct_positions_by_overlap(target_features, context_features, overlap, maxoverlap, seqlength, seqname):
	#overlap, maxoverlap, seqlength = [args.overlap, args.overlap_maxdist, len(seq_record.seq)]
	
	context_overdist = set()
	
	# Check context_features all match in positions
	for name, feats in context_features.items():
		locations = [feat.location for feat in feats]
		if(not all(locations[0] == loc for loc in locations)):
			err = "Warning, positions of " + str(len(locations)) + " annotations for " + name + " in " + seqname + " do not match. If these are multiple distinct loci, this should be fine, but if these should cover the same locus, you may get incorrect results.\n"
			for i, feat in enumerate(feats):
				err += "\t(" + str(i+1) + ") " + feat.type +" is located at bases " + str(int(feat.location.start)+1) + " to " + str(int(feat.location.end)) + "\n"
			sys.stderr.write(err)
	
	# Work through combinations of target and context features
	
	for target in target_features:
		#target = target_features[0]
		tpos = [int(target.location.start), int(target.location.end)]
		for context_name in context_features:
			#context_name = list(context_features.keys())[0]
			context = context_features[context_name][0]
			cpos = [int(context.location.start), int(context.location.end)]
			
			#tpos, cpos,overlap = [[1,6],[3,4],-1]
			
			# Set orientation (+ve, context follows target)
			orientation = 0
			# Set current distance between selected positions
			distance = 0
			# Set the index of the target position to change
			target_tpos_i = None
			
			# Find the structure of the overlap
			if(min(cpos)-min(tpos) == max(tpos)-max(cpos)):
				# Current positions are completely even, can't figure it out!
				print("error")
				continue
			elif(min(tpos) < min(cpos)):
				# Target is before context, overlap should be latter position of target and first position of context
				orientation = 1
				distance = min(cpos)-max(tpos)
				target_tpos_i = tpos.index(max(tpos))
			elif(max(tpos) > max(cpos)):
				# Target is after context, overlap should be first position of target and latter position of context
				orientation = -1
				distance = max(cpos)-min(tpos)
				target_tpos_i = tpos.index(min(tpos))
			else:
				# Target is completely within context, can't figure it out!
				print("error")
				continue
			
			if(abs(distance) > maxoverlap + overlap[context_name]):
				context_overdist.add(context_name)
				continue
			
			# Calculate the exact new position
			corrected_tpos = tpos[target_tpos_i] + distance + orientation * (overlap[context_name])
			corrected_tpos
			# Ensure the position is not outside the contig
			corrected_tpos = 0 if corrected_tpos < 0 else corrected_tpos
			corrected_tpos = seqlength if corrected_tpos > seqlength else corrected_tpos
			
			# Overwrite the relevant target end position
			if(target_tpos_i == 0):
				target.location = SeqFeature.FeatureLocation(corrected_tpos, target.location.end, target.location.strand)
			else:
				target.location = SeqFeature.FeatureLocation(target.location.start, corrected_tpos, target.location.strand)
	return(context_overdist)

def correct_feature_by_query(feat, query_spec, seq_record, seqname, distance, featurename):
	#query_spec, distance, featurename = [stringspec, args.search_distance, args.annotation[0]]
	
	feat_start, feat_finish = feat.location.start, feat.location.end
	errstart = "Warning: sequence " + seqname + " has "
	
	for end, search in query_spec.items():
		#end, search = next(iter(query_spec.items()))
		
		# Unpack search tuple
		code, query, out_rf, selector = search if len(search) == 4 else search + tuple("X")
		selector = out_rf if code == 'A' else selector
		
		# Check if already ends with the searched sequence - REMOVED AS EXISTING SEQUENCE MAY BE SHORTER SUBSET OF DESIRED SEQUENCE e.g. TA TAA
		#if(end_already_correct(feat.extract(seq_record.seq), query, end, code, out_rf)):
		#	continue
		
		# Generate sequence for searching
		subject_sequence, subject_start = extract_subject_region(seq_record, feat, end, code, distance)
		
		# Find locations of query
		results = dict()
		for q in query.split("/"):
			# Work through locations of hits
			for i in list(find_all(str(subject_sequence), q)):
				# If the current hit location exists and is longer than the current hit, do not change, else add current hit length
				results[i] = results[i] if i in results.keys() and len(q) < results[i] else len(q)
		
		# Remove results if selector is N, FC or LC
		errmid = ""
		if(len(results) == 0):
			errmid = "no matches of "
		if(selector == "N" and len(results) > 1):
			errmid = "multiple matches of "
			results = {}
		if(selector == "FC"):
			results = {i:l for i, l in results.items() if i <= distance}
			errmid = "no first closest matches of "
		elif(selector == "LC"):
			results = {i:l for i, l in results.items() if i >= distance}
			errmid = "no last closest matches of "
		
		# Retain only locations in the specified reading frame
		if(code == 'N' and out_rf != "*" and len(results) > 0):
			if(end == "start"):
				results = {i:l for i, l in results.items() if (i + distance) % 3 + 1 == int(out_rf)}
				# Retain location if that location's rf (l+1)%3 is equal to the (target rf converted to subject rf)
			else:
				results = {i:l for i, l in results.items() if (abs(feat.location.end - feat.location.start) - distance + i - 1) % 3 + 1 == int(search[2])}
			errmid = "no matches in the specified frame of " if len(results) == 0 else errmid
		
		truncated = False
		codon_start = 1 if end == "start" else None
		
		if(end == "start"):
			
			# First check if truncated
			start_distance = len(seq_record) - feat.location.end if feat.location.strand == -1 else feat.location.start
			truncated = start_distance < distance
			
			# Retain only start locations that generate realistic amino acid sequences
			results = {i:l for i, l in results.items() if is_inframe(i, l, feat.location.strand, code, end, distance, subject_start, feat_start, feat_finish, seq_record, args.translation_table)}
			
			# If no feasible results at the start position, instead find the closest in-frame position 
			if(len(results) == 0):
				errmid = "no ORF-producing matches (will set to closest ORF) of "
				
				# Set the current position - if normal, this is the current start position, if truncated, this is the end of the contig
				contig_start = list(find_all(str(subject_sequence), 'N'))[-1] + 1 if truncated else None
				current_position = contig_start if truncated else distance
				
				# Set the correction - if normal, this is 0, if truncated, this is 1, to ensure no results outside the contig
				correction = 1 if truncated else 0
				
				# Set up the three alternative results
				results = { current_position + v + correction : 1 for v in [-1, 0, 1] }
				
				# Find the result with suitable ORF
				results = {i:l for i, l in results.items() if is_inframe(i, l, feat.location.strand, code, end, distance, subject_start, feat_start, feat_finish, seq_record, args.translation_table)}
				
				# If truncated, set result to contig start but note codon position
				if(len(results) > 0 and truncated):
						codon_start = sorted(results.keys())[0] - contig_start + 1
						results = {contig_start : 1}
				
			
		else:
			# Prioritise longer matches by removing any matches shorter than the longest match
			max_length = max(results.values())
			results = {i:l for i, l in results.items() if l == max_length}
			
			# If searching for finish string and the annotation is likely truncated, remove any incomplete stop codons
			
			# Find the distance to the finish of the contig from the current position (strand-dependent)
			finish_distance = feat.location.start if feat.location.strand == -1 else len(seq_record) - feat.location.end
			
			# Remove partial stops if annotation likely truncated
			if(finish_distance < distance):
				truncated = True
				results = {i:l for i,l in results.items() if l > 2}
				# If no results remain, set result to be the end of the contig
				if(len(results) == 0):
					errmid = "no >2 matches near truncated finish (will set to contig end) of "
					results = {str(subject_sequence).find('N') - 1 : 1}
		
		# Parse location results
		errend = query  + " at the " + end + " of " + featurename + "\n"
		
		if(len(results) > 0):
			
			# Select the first location by default (which also is the LC location if LC set)
			locations = sorted(results.keys())
			location = locations[0]
			
			if(len(locations) > 1):
				if(selector in ['L', 'FC']):
					# If more than 1 locations but the user wants the last or first closest selected
					location = locations[-1]
				elif(selector == 'C'):
					loc_dist = [abs(distance-l) for l in locations]
					if(loc_dist.count(min(loc_dist)) == 1):
						location = locations[loc_dist.index(min(loc_dist))]
					else:
						errmid = "multiple closest matches (taking first) of "
			
			feat_start, feat_finish = get_newends(location, results[location], feat.location.strand, end, distance, code, subject_start, feat_start, feat_finish, truncated)
		else:
			errmid = "no succesful matches of " if errmid == "" else errmid
		
		if(errmid != ""):
			sys.stderr.write(errstart + errmid + errend)
	
	return(feat_start, feat_finish, codon_start)

def get_newends(location, length, strand, end, distance, code, subject_start, feat_start, feat_finish, truncated):
	
	# Convert location if on reverse strand
	location = location if strand == 1 else abs(location - (2*distance + 1))
	
	# Generate the new end position
		# Correct by length of the match if at the finish end
	change = location + strand * length if end == "finish" else location
		# Multiply by 3 if AA
	change = change * 3 if code == 'A' else change
		# Calculate
	newend = subject_start + change
	
	# Apply new end to appropriate end
	if((end == "start" and strand == 1) or (
			end == "finish" and strand == -1)):
		feat_start = SeqFeature.BeforePosition(newend) if truncated else SeqFeature.ExactPosition(newend)
	else:
		feat_finish = SeqFeature.AfterPosition(newend) if truncated else SeqFeature.ExactPosition(newend)
	
	return(feat_start, feat_finish)

def stopcount(seq_record, table, frame = (1,2,3)):
	
	# Check input types
	run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
	
	# Run counting
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', BiopythonWarning)
		counts = [seq_record.seq[(i-1):].translate(table = table).count("*") for i in run_frame]
	
	# Return string or list depending on length
	if(len(counts) > 1):
		return counts
	else:
		return counts[0]

def is_inframe(location, length, strand, code, end, distance, subject_start, feat_start, feat_finish, seq_record, table):
	#location, length, strand, table = [19, results[20], feat.location.strand, args.translation_table]
	# Generate the potential end position for this location
	potstart, potfinish = get_newends(location, length, strand, end, distance, code, subject_start, feat_start, feat_finish, False)
	
	# Build the potential new feature for this location
	potfeat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(potstart, potfinish), strand = strand)
	
	# Extract the sequence for this potential new feature
	potseq = SeqRecord.SeqRecord(potfeat.extract(seq_record.seq))
	
	# Count the stops in this sequence and return true if less than or equal to 1
	return(stopcount(potseq, table, 1) <= 1)



def extract_subject_region(seqrecord, feat, end, code, distance):
	'''For a given feature and end ("start" or "finish"), extract a sequence of either nucleotides (code = 'N') or amino acids (code = 'A'), consisting of the first or last position plus or minus positions equal to distance in the reading direction of the feature'''
	#seqrecord = seq_record
	
	# Find the centre point and distances for the subject region
	featend = "featstart"
	central_position = feat.location.start
	distances = (distance, distance + 1)
	if((end == "start" and feat.location.strand == -1) or (end == "finish" and feat.location.strand == 1)):
		featend = "featfinish"
		central_position = feat.location.end
		distances = (distance + 1, distance)
	
	# Convert for amino acids
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
	
	# Delimit the region
	end_positions = (central_position - distances[0], central_position + distances[1])
	
	# Truncate the region if it exceeds the contig
	start_position = 0 if end_positions[0] < 0 else end_positions[0]
	finish_position = len(seqrecord) if end_positions[1] > len(seqrecord) else end_positions[1]
	
	# Truncate the region if it exceeds the other end of the annotation
	if(featend == "featstart"):
		finish_position = feat.location.end if finish_position > feat.location.end else finish_position
	else:
		start_position = feat.location.start if start_position < feat.location.start else start_position
	
	# Generate the feature	
	subject_feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start_position, finish_position), strand = feat.location.strand)

	# Extract the sequence
	sequence = subject_feat.extract(seqrecord.seq)
	correction = [start_position - end_positions[0], end_positions[1] - finish_position]
	correction = correction[::-1] if feat.location.strand == -1 else correction
	sequence = correction[0] * "N" + sequence + correction[1] * "N"
	
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
		return(any(end == "start" and frame == "1" and nuc_seq.startswith(q) or # Starts with sequence, rf is 1
			   (end == "finish" and nuc_seq.endswith(q) and nuc_seq.rfind(q) % 3 + 1 == int(frame) )) for q in query_seq.split("/"))
	else:
		aa_seq = nuc_seq.translate(table = args.translation_table)
		return(any((end == "start" and aa_seq.startswith(q)) or
		           (end == "finish" and aa_seq.endswith(query_seq)) for q in query_seq.split("/")))

if __name__ == "__main__":
	
	# Read in arguments
	
	#args = parser.parse_args(['-a', "ATP6", '-c', 'ATP8', '-o', '7', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-04_2edited/SPSO00168.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "CYTB", '-c', 'TRNS(UGA)', '-o', '2', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00109.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "NAD2", '-o', 'TRNW,2', '-o', 'TRNS,-20', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00001.gb', '-m', '50'])
	#args = parser.parse_args(['-a', "NAD", '-c', 'TRNH(GUG)', '-o', '0', '-i', '/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00409.gb', '-m', '50'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/QINL005.gb', '-a', 'ND5', '-s', 'N,ATT/ATA/ATG/ATC,1,L', '-d', '3', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/testing/BIOD00550.gb', '-a', 'ND2', '-s', 'N,ATA/ATG/ATC/TTG/ATT,*,LC', '-d', '20', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/CCCP00094.gb', '-a', 'NAD1', '-f', 'N,TAA/TAG,1,F', '-d', '220', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-04_2edited/BIOD00622.gb', '-a', 'ATP6', '-f', 'N,TAG/TAA/TA,1,F', '-d', '15', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00010.gb', '-a', 'NAD4', '-s', 'N,ATG/ATA,1,F', '-d', '6', '-t', '5'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/CCCP00017.gb', '-y', 'gene'])
	#args = parser.parse_args(['-i','/home/thomas/Documents/NHM_postdoc/MMGdatabase/gbmaster_2020-02-12_2edited/BIOD00109.gb', '-a', 'ND6', '-s', 'N,ATT/ATA/ATC/TTG/TTT,*,FC', '-d', '6', '-t', '5'])
	
	args = parser.parse_args()
	
	# Check arguments
	
	stringspec = dict()
	overlap = {}
	if(args.overlap is not None):
		if(len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --overlap")
		if(args.startstring or args.finishstring):
			sys.exit("Error: --startstring or --finishstring not compatible with --overlap. Run consecutive iterations to perform both options")
		sys.stderr.write("Running overlap autocorrection for %s based on:\n"  % (args.annotation[0]))
		for overlapstring in args.overlap:
			#overlapstring = args.overlap[0]
			overlap_list = overlapstring.split(",")
			err = "Error: overlap string " + overlapstring
			if(len(overlap_list) != 2):
				sys.exit(err + " does not contain exactly two comma-separated elements")
			elif(not str_is_int(overlap_list[1])):
				sys.exit(err + " does not have an integer as the second item")
			sys.stderr.write("\toverlap of %s bp with %s\n" % tuple(reversed(overlap_list)))
			overlap[overlap_list[0]] = int(overlap_list[1])
	elif(args.startstring is not None or args.finishstring is not None):
		if(args.annotation is None or len(args.annotation) != 1):
			sys.exit("Error: please supply one and only one annotation name to --annotation if using --*string")
		
		sys.stderr.write("Running string searching for %s using:\n" % (args.annotation[0]))
		for end, searchstring in zip(['start','finish'], [args.startstring, args.finishstring]):
			if(searchstring is None):
				continue
			if(end is 'start' and args.translation_table is None):
				sys.exit("Error: --translation_table must be specified if searching for --startstring")
			sys.stderr.write("\t%s: %s\n" % (end, searchstring))
			search = tuple(searchstring.split(","))
			err = "Error: search string " + searchstring
			if(search[0] not in ['A', 'N']):
				sys.exit(err + " has an unrecognised sequence type")
			elif(search[0] == 'A'):
				if(len(search) != 3):
					sys.exit(err + " does not have only two or three items")
				elif(any(s not in list("/GPAVLIMCFYWHKRQNEDST*") for s in search[1])):
					sys.exit(err + " is specified as amino acid but non-standard character(s) included (must be GPAVLIMCFYWHKRQNEDST*)")
				elif(str_is_int(search[2])):
					sys.exit(err + " is searching for an amino acid sequence but includes what seems to be a reading frame")
				elif(args.translation_table is None):
					sys.exit(err + " is specified as amino acid but no --translation_table given")
			else:
				if(len(search) != 4):
					sys.exit(err + " does not have only three or four items")
				elif(any(s not in list("/ATGC") for s in search[1])):
					sys.exit(err + " is specified as nucleotide but non-standard character(s) included (must be ATGC)")
				elif(search[2] not in ['1','2','3','*']):
					sys.exit(err + " has an unrecognised reading frame")
			
			mmhi = search[3] if search[0] == 'N' else search[2]
			if( mmhi not in ['F','L','C','FC','LC','N'] ):
				sys.exit(err + " has an unrecognised multiple-hit handling instruction ( must be F, FC, C, LC, L or N)")
			
			stringspec[end] = search
		
	elif(args.syncronise is not None):
		if(args.annotation is not None or args.overlap is not None):
			sys.exit("Error: no --annotation or --overlap should be supplied to --syncronise")
		if(args.syncronise not in ['gene', 'CDS', 'tRNA']):
			sys.exit("Erorr: value passed to --syncronise should be gene, CDS or tRNA")
		sys.stdout.write("Running position syncronisation on %s annotations\n" % (args.syncronise))
	else:
		sys.exit("Error: please supply a value to --overlap, to --startstring and/or --endstring, or to --syncronise")
	
	# Read and parse gene name variants
	
	namevariants = loadnamevariants()
	
	# Find universal names for inputs
	err = "Error: unrecognised locus name supplied to"
	
	if(args.syncronise is not None):
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
		
		# Get features and parse errors
		feature_names = args.annotation
		feature_names += overlap.keys() if overlap is not None else []
		
		features, record_unrecognised_names, record_unidentifiable_features = get_features_from_names(seq_record, feature_names, namevariants)
		unrecognised_names.update(record_unrecognised_names)
		
		if(len(record_unidentifiable_features) > 0):
			unidentifiable_features[seqname] = record_unidentifiable_features
		
		if(args.overlap is not None):
			target_features = features.pop(args.annotation[0]) if args.annotation[0] in features.keys() else []
			context_features = features
			
			ntf = len(target_features)
			ncf = sum([len(cfl) for gene, cfl in context_features.items()])
			if(ntf > 0 and ncf > 0):
				record_context_overdist = correct_positions_by_overlap(target_features, context_features, overlap, args.overlap_maxdist, len(seq_record.seq), seqname)
				if(len(record_context_overdist) > 0):
					context_overdist[seqname] = record_context_overdist
			else:
				if(ntf == 0): missing_annotation.add(seqname)
				if(ncf == 0): missing_context.add(seqname)
				#err = "Warning, sequence " + seqname + " has"
				#if ntf == 0: err += " no " + args.annotation[0] + " annotation" 
				#if ntf == 0 and ncf == 0: err += " and"
				#if ncf == 0: err += " none of the specified context annotation(s)"
				#sys.stderr.write(err + "\n")
		
		elif(len(stringspec) > 0):
			features = features[args.annotation[0]]
			
			nf = len(features)
			
			if(nf > 0):
				for feat in features:
					#feat = features[0]
					#feat.extract(seq_record.seq)
					
					# Get new start and finish positions
					corrected_start, corrected_finish, codon_start = correct_feature_by_query(feat, stringspec, seq_record, seqname, args.search_distance, args.annotation[0])
					
					#sys.stderr.write("Before: %i , %i; After: %i, %i\n" % (int(feat.location.start), int(feat.location.end), corrected_start, corrected_finish))
					
					if(corrected_start == feat.location.start and corrected_finish == feat.location.end):
						continue
					else:
						feat.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, feat.location.strand)
						if(codon_start is not None):
							feat.qualifiers['codon_start'] = codon_start
			else:
				#err = "Warning, sequence " + seqname + " has no " + args.annotation[0] + " annotation(s)\n"
				#sys.stderr.write(err)
				missing_annotation.add(seqname)
		
		elif(args.syncronise is not None):
			syncronise_features(features, args.syncronise, seqname)
		
		output_records.append(seq_record) 
	
	if len(output_records)>0 : SeqIO.write(output_records, sys.stdout, "genbank")
	
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
