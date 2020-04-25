#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 09:46:36 2020

@author: thomas
"""

import sys
import copy
import urllib.request
import re
from statistics import mode #, stdev

from collections import defaultdict, Counter
from Bio import SeqFeature, SeqRecord, AlignIO
from Bio.Align import AlignInfo

from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
	warnings.simplefilter('ignore', BiopythonWarning)




def cumsum(lis):
	total = 0
	for i in lis:
		total += i
		yield total

def ungapped_distance(seq_str, locations, end):
	# Reverse the sequence if we're looking at the finish
	sequence = seq_str[::-1] if end is 'finish' else seq_str
	
	# Extract the in between sequence
	between_sequence = sequence[min(locations):max(locations)]
	
	# Find the ungapped distance	
	return(len(between_sequence) - between_sequence.count("-"))

def ishigher(x, y):
	for v in y:
		if(v > x):
			return(False)
	return(True)

def find_maxima(lis, d):
	#lis = seg
	#d = 10
	maxima = []
	i = 0
	while i < len(lis):
		s, f = [i - d, i + d]
		s = s if s >= 0 else 0
		if(ishigher(lis[i], lis[s:f])):
			maxima.append(i)
			i += d + 1
		else:
			i += 1
	return(maxima)

def find_0intercepts(lis):
	intercepts = set()
	i=0
	while i < len(lis):
		pair = lis[i:i+2]
		if(min(pair) <= 0 and max(pair) >= 0):
			abspair = [abs(p) for p in pair]
			intercepts.add(i + abspair.index(min(abspair)))
		i += 1
	return(sorted(list(intercepts)))

def parse_alignment(path, frame):
	#path, frame = [args.match_alignment, args.force_alignment_frame]
	frame = frame - 1 if frame is not None else None
	
	alignment = AlignIO.read(path, "fasta")
	
	# Get modal start and finish positions
	prelim = dict()
	for seq_record in alignment:
		prelim[seq_record.id] = {e: len(re.search(r, str(seq_record.seq)).group()) for e, r in zip(["start", "finish"],["^-*", "-*$"])}
	start, finish = zip(*[[d[e] for e in ["start", "finish"]] for d in prelim.values()])
	
	modes = dict()
	for e, v in zip(['start', 'finish'], [start, finish]):
		modes[e] = mode(v) if frame is None else round(mode(v)/3)*3+frame
	
	# Generate consensus
	alignment_summary = AlignInfo.SummaryInfo(alignment)
	consensus = str(alignment_summary.gap_consensus())
	#modal_consensus = consensus[modes['start']:-modes['finish']]
	
	# Build consensus gap pattern and find positions for alignment body
	# TERMINOLOGY:
	# ----X--X---XXXXXXXXXXXX--X-X-----
	#-----X--X---XXXXXXXXXXXX--X-X--X--
	#-X---X--X---XXXXXXXXXXXX--X-X-----
	#     1      2          3    4     
	# 1: modal start
	# 2: body start
	# 3: body finish
	# 4: modal finish
		
		# Find ratio of characters to gaps
	def gapratio(seqstring):
		ngaps = seqstring.count('-')
		nchars = len(seqstring) - ngaps
		return(round(ngaps/nchars,5))
	
	consensus_gapratio = gapratio(consensus)
	#modal_gapratio = gapratio(modal_consensus)
	
		# Find pattern of gaps
	#modal_gap_pattern = [1 if b is '-' else -modal_gapratio for b in modal_consensus]
	consensus_gap_pattern = [1 if b is '-' else -consensus_gapratio for b in consensus]
	
	c_modes = dict()
	for e, m in modes.items():
		#e, m = list(modes.items())[0]
		
		# Sort the gap pattern in the direction required.
		gp = reversed(consensus_gap_pattern) if e == 'finish' else consensus_gap_pattern
		
		# Calculate the cumulative sum (add 0 to end in case rounding errors prevent return to 0), then take the segment between the modal position and the next time the ratio of gaps:characters passes 0
		cusum = list(cumsum(gp)) + [0]
		
		#from matplotlib import pyplot
		#pyplot.plot(cusum)
		#print(cusum)
		
		seg = cusum[m:]
		body, dist = [0, 0]
		if (m > 0 and any([s > 0 for s in seg]) and find_0intercepts(seg)[0] > 0):
			segshort = seg[:find_0intercepts(seg)[0]]
			#pyplot.plot(seg)
			#pyplot.plot(segshort)
			
			# Find the location of the start of the body of the alignment with multiple methods
			body = []
				# Find the maximum gap position
			body.append(segshort.index(max(segshort)))
				# Find the first regional maximum of the gap pattern; 
			regional_max = [x for x in find_maxima(segshort, 30) if x != 0]
			body.append(regional_max[0] if len(regional_max)>0 else None)
			
				# Pick the first
			body = min([b for b in body if b is not None])
			
			# add this to the modal position and find the closest amino acid location
			body = round((m+body)/3)*3+frame
			
			# Calculate the number of ungapped consensus positions between the mode and this position
			dist = ungapped_distance(consensus, [modes[e], body], e)
			
		# Output these values
		c_modes[e] = (body, -dist)
	
	modes = {e:(v, 0) for e, v in modes.items()}
	
	# Find approximate ungapped distance to target position for each sequence
	def get_seqdists(alignment, targets):
		output = dict()
		
		for seq_record in alignment:
			#seq_record = [sr for sr in alignment if sr.id == 'BIOD02328'][0]
			dists = dict()
			sequence = str(seq_record.seq)
			
			for e in ['start', 'finish']:
				#e = 'start'
				#e = 'finish'
				# Set up the two locations in a list
				locations = [prelim[seq_record.id][e], targets[e][0]]
				
				#Check whether the sequence location is inside the target
				check_seq = sequence
				correction = targets[e][1]
				if(locations[0] > locations[1]):
					check_seq = consensus
					correction = -correction
				
				# Get the ungapped distance modified by the standard distance
				dists[e] = ungapped_distance(check_seq, locations, e) + correction
				
				# Correct the distance if inside the body
				dists[e] = dists[e] if(locations[0] < locations[1]) else dists[e] * -1
				
				# Correct the distance depending on end
				dists[e] = dists[e] * -1 if e == 'finish' else dists[e]
				
			output[seq_record.id] = dists
		
		return(output)
	
	
	output = {'mode': get_seqdists(alignment, modes),
			 'body': get_seqdists(alignment, c_modes)}
	output = {s:{m:output[m][s] for m in output.keys()} for s in output['mode'].keys()}
	
	# Find standard deviation of the distance
	#start, finish = zip(*[[d[e] for e in ["start", "finish"]] for d in output.values()])
	
	#std_dev = {e:round(stdev(v)) for e, v in zip(["start", "finish"], [start, finish])}
	
	return(output)


def check_arguments(arguments):
	stringspecdict = dict()
	overlapset = {}
	alignment_dists = dict()
	alignment_sds = []
	if(arguments.syncronise is not None):
		if(arguments.startstring is not None or
			arguments.finishstring is not None or
			arguments.overlap is not None or 
			arguments.match_alignment is not None):
			sys.exit("Error: --syncronise is not compatible with --startstring, --finishstring, --annotation, --overlap or --match_alignment")
		
		if(arguments.syncronise not in ['gene', 'CDS', 'tRNA']):
			sys.exit("Error: value passed to --syncronise should be gene, CDS or tRNA")
			sys.stdout.write("Running position syncronisation on %s annotations\n" % (arguments.syncronise))
		
	elif(arguments.annotation is not None):
		
		if(arguments.overlap is not None):
			
			if(arguments.match_alignment is None):
				overlapset = parse_overlap(arguments.overlap, arguments.annotation)
			else:
				sys.exit("Error: --match_alignment and --overlap are mutually exclusive")
			
		elif(arguments.match_alignment is not None):
			
			# Load in and process the alignment if doing matching
			alignment_dists = parse_alignment(arguments.match_alignment, arguments.force_alignment_frame)
			sys.stderr.write("Completed loading and parsing "+ arguments.match_alignment + "\n")
			
		elif(arguments.startstring is None and arguments.finishstring is None):
			
			sys.exit("Error: insufficient arguments - are you missing at least one of --overlap, --match_alignment, --startstring and/or --finishstring?")
		
		
		if(arguments.startstring is not None or arguments.finishstring is not None):
			stringspecdict = parse_stringsearch(arguments.annotation, arguments.startstring, arguments.finishstring, arguments.translation_table, arguments.match_alignment is not None)
		
	else:
		sys.exit("Error: insufficient arguments - are you missing --annotation?")
	
	return(stringspecdict, overlapset, alignment_dists)

def standardise_names(arguments, overlapdict, namevariants):
	err = "Error: unrecognised locus name supplied to"
	
	if(arguments.syncronise is not None):
		arguments.annotation = list(set(namevariants.values()))
	else:
		if(arguments.annotation is not None):
			if(arguments.annotation.upper() in namevariants):
				arguments.annotation = namevariants[arguments.annotation.upper()]
			else:
				err = err + " --annotation"
				sys.exit(err)
		
		if(overlapdict is not None):
			if all(c.upper() in namevariants for c in overlapdict.keys()):
				overlapdict = { namevariants[c.upper()]:o for c,o in overlapdict.items() }
			else:
				#die with error
				err = err + " --overlap"
				sys.exit(err)
	
	return(arguments, overlapdict)



def str_is_int(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False

def parse_overlap(overlap_in, annotation_in):
	sys.stderr.write("Running overlap autocorrection for %s based on:\n"  % (annotation_in))
	output = {}
	
	for overlapstring in overlap_in:
		#overlapstring = args.overlap[0]
		overlap_list = overlapstring.split(",")
		err = "Error: overlap string " + overlapstring
		if(len(overlap_list) != 2):
			sys.exit(err + " does not contain exactly two comma-separated elements")
		elif(not str_is_int(overlap_list[1])):
			sys.exit(err + " does not have an integer as the second item")
		sys.stderr.write("\toverlap of %s bp with %s\n" % tuple(reversed(overlap_list)))
		output[overlap_list[0]] = int(overlap_list[1])
	
	return(output)

def parse_stringsearch(annotation_in, startstring_in, finishstring_in, translation_table, match_alignment):
	#annotation_in, startstring_in, finishstring_in, translation_table, match_alignment = [args.annotation[0], args.startstring, args.finishstring, args.translation_table, args.match_alignment is not None]
	output = dict()
	sys.stderr.write("Running string searching for %s using:\n" % (annotation_in))
	
	for end, searchstring in zip(['start','finish'], [startstring_in, finishstring_in]):
		if(searchstring is None):
			continue
		if(end is 'start' and translation_table is None):
			sys.exit("Error: --translation_table must be specified if searching for --startstring")
		sys.stderr.write("\t%s: %s\n" % (end, searchstring))
		
		search = searchstring.split(",")
		
		err = "Error: search string " + searchstring
		if(search[0] not in ['A', 'N']):
			sys.exit(err + " has an unrecognised sequence type")
		elif(search[0] == 'A'):
			if((match_alignment and len(search) !=2) or (not match_alignment and len(search) != 3)):
				sys.exit(err + " does not have the correct number of items")
			elif(any(s not in list("/GPAVLIMCFYWHKRQNEDST*") for s in search[1])):
				sys.exit(err + " is specified as amino acid but non-standard character(s) included (must be GPAVLIMCFYWHKRQNEDST*)")
			elif(not match_alignment and str_is_int(search[2])):
				sys.exit(err + " is searching for an amino acid sequence but includes what seems to be a reading frame")
			elif(translation_table is None):
				sys.exit(err + " is specified as amino acid but no --translation_table given")
		else:
			if((match_alignment and len(search) !=3) or (not match_alignment and len(search) != 4)):
				sys.exit(err + " does not have the correct number of items")
			elif(any(s not in list("/ATGC") for s in search[1])):
				sys.exit(err + " is specified as nucleotide but non-standard character(s) included (must be ATGC)")
			elif(search[2] not in ['1','2','3','*']):
				sys.exit(err + " has an unrecognised reading frame")
		
		if(not match_alignment):
			mmhi = search[3] if search[0] == 'N' else search[2]
			if( mmhi not in ['F','L','C','FC','LC','N'] ):
				sys.exit(err + " has an unrecognised multiple-hit handling instruction ( must be F, FC, C, LC, L or N)")
		else:
			search = search + ['C']
		
		output[end] = search
	
	return(output)

def find_all(a_str, sub):
	start = 0
	while True:
		start = a_str.find(sub, start)
		if start == -1: return
		yield start
		start += 1

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
	
	names = [names] if isinstance(names, str) else names
	
	features = defaultdict(list)
	unrecognised_names = set()
	unidentifiable_features = set()
	
	for feat in seqrecord.features:

		# Remove any translations
		if('translation' in feat.qualifiers.keys()):
			del(feat.qualifiers['translation'])
		
		# Extract the tag that contains the feature name
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
			continue
		
		# Find the standard name
		if(featname in namevariants):
			name = namevariants[featname]
			if(name in names):
				features[name].append(feat)
		else:
			unrecognised_names.add(featname)
	
	return(features, unrecognised_names, unidentifiable_features)

def syncronise_features(name, feats, synctype, seqname):
	#synctype = 'gene'
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
		return
	if(len(other_feats) < 1):
		#sys.stderr.write(warnstart + "non-" + synctype + warnend)
		return
	
	# Check that the other features are all the same
	if(not all(other_feats[0].location == feat.location for feat in other_feats)):
		sys.stderr.write("Warning, positions of " + str(len(other_feats)) + " non-" + synctype + " annotations for " + name + " in " + seqname + " do not match, this entry will not be modified\n")
		return
		
	# Correct the target features
	for feat in target_feats:
		feat.location = other_feats[0].location

def correct_positions_by_overlap(target, context_features, overlap, maxoverlap, seqlength, seqname):
	#target, maxoverlap, seqlength = [feat, args.maxdist, len(seq_record.seq)]
	
	context_overdist = set()
	outfeat = copy.deepcopy(target)
	corrected_start, corrected_finish = [target.location.start, target.location.end]
	
	# Work through combinations of target and context features
	
	tpos = [int(target.location.start), int(target.location.end)]
	for context_name, context_feats in context_features.items():
		#context_name, context_feats = list(context_features.items())[0]
		for context in context_feats:
			#context = context_feats[1]
			cpos = [int(context.location.start), int(context.location.end)]
			
			# Find gap between context and target
			gap = [min(cpos)-max(tpos), min(tpos)-max(cpos)]
			gap = [g for g in gap if g > 0]
			
			# Check gap is permissible
			if(len(gap) > 0 and gap[0] > maxoverlap + overlap[context_name]):
				context_overdist.add(context_name)
				continue
			
			# Set orientation (+ve, context follows target)
			orientation = 0
			# Set current distance between selected positions
			distance = 0
			# Set the index of the target position to change
			target_tpos_i = None
			
			
			# Find cross-wise distances
			crossdist = [max(cpos) - min(tpos), min(cpos) - max(tpos)]
			
			# Find the structure of the overlap
			if(abs(crossdist[0]) == abs(crossdist[1])):
				sys.stderr.write("Warning: orientation of " + context_name + " and " + "target annotations completely match, cannot perform overlap correction\n")
				continue
			elif(abs(crossdist[1]) < abs(crossdist[0])):
				# Target is before context, overlap should be latter position of target and first position of context
				orientation = 1
				distance = crossdist[1]
				target_tpos_i = tpos.index(max(tpos))
			else:
				# Target is after context, overlap should be first position of target and latter position of context
				orientation = -1
				distance = crossdist[0]
				target_tpos_i = tpos.index(min(tpos))
			
			# Calculate the exact new position
			corrected_tpos = tpos[target_tpos_i] + distance + orientation * (overlap[context_name])
			
			# Ensure the position is not outside the contig
			corrected_tpos = 0 if corrected_tpos < 0 else corrected_tpos
			corrected_tpos = seqlength if corrected_tpos > seqlength else corrected_tpos
			
			# Check that the orientation hasn't been flipped
			if((target_tpos_i == 0 and corrected_tpos >= int(target.location.end)) or (target_tpos_i == 1 and corrected_tpos <= int(target.location.start))):
				continue
			
			# Overwrite the relevant target end position
			if(target_tpos_i == 0):
				corrected_start = corrected_tpos
			else:
				corrected_finish = corrected_tpos
	
	outfeat.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, outfeat.location.strand)
	
	return(outfeat, context_overdist)

def correct_location_if_valid(location, correction, checkvalue):
	if(str_is_int(str(location)) and location != checkvalue):
		return(location + correction)
	else:
		return(location)

def correct_feature_by_alignment(feat, query_spec, distances, featname, seqname, seqlength):
	#query_spec, distances, featname, seqlength = [stringspec, alignment_distances[seqname], name, len(seq_record)]
	
	outfeat = copy.deepcopy(feat)
	
	for end in ['start', 'finish']:
		#end = 'start'
		#end = 'finish'
		
		# Skip if not processing this end or if the end is already correct
		if(end not in query_spec.keys() or (end == 'start' and distances[end] == 0)):
			continue
		
		locations = [outfeat.location.start, outfeat.location.end]
		dist = distances[end]
		# Set finish distance to the first base of the stop to prevent double stops
		dist -= [2,0,1][(locations[1]-locations[0])%3] if end == 'finish' else 0
		
		if((end == 'start' and feat.location.strand == 1) or 
		   (end == "finish" and feat.location.strand == -1)):
			
			locations[0] =  correct_location_if_valid(locations[0], feat.location.strand * dist, 0)
			
		else:
			
			locations[1] =  correct_location_if_valid(locations[1], feat.location.strand * dist, seqlength)
		
		if( locations[0] > locations[1] ):
			sys.stderr.write("Warning: changing " + end + " of " + featname + " in " + seqname + " to match alignment causes incorrect orientation, no change made\n")
		elif( locations[0] < 0 or locations[1] > seqlength):
			sys.stderr.write("Warning: changing " + end + " of " + featname + " in " + seqname + " to match alignment causes annotation to exceed contig, no change made\n")
		else:
			outfeat.location = SeqFeature.FeatureLocation(locations[0], locations[1], outfeat.location.strand)
		
	movement = [outfeat.location.start - feat.location.start, outfeat.location.end - feat.location.end]
	movement = [-m for m in reversed(movement)] if feat.location.strand == -1 else movement
	distances_moved = {e:v for e, v in zip(['start', 'finish'], movement)}
	
	return(outfeat, distances_moved)

def get_newends(location, length, feat_start, feat_finish, strand, end, distance, code, subject_start, truncated):
	
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

def stopcount(seqr, table, frame = (1,2,3), includefinal = True):
	# seqr, frame, includefinal = [potseq, 1, False]
	
	# Check input types
	run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
	
	# Run counting
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', BiopythonWarning)
		aa = [str(seqr.seq[(i-1):].translate(table = table)) for i in run_frame]
		aa = [re.sub("\*+$", "", a) for a in aa] if not includefinal else aa
		counts = [a.count("*") if len(a) > 0 else 10 for a in aa]
	
	# Return string or list depending on length
	if(len(counts) > 1):
		return counts
	else:
		return counts[0]

def get_stopcounts(location, length, feat, code, end, distance, subject_start, seq_record, table):
	#location, length, strand, table = [19, results[20], feat.location.strand, args.translation_table]
	#location, length, strand, table = list(results.items())[1]+(feat.location.strand, translation_table)
	# Generate the potential end position for this location
	potstart, potfinish = get_newends(location, length, feat.location.start, feat.location.end, feat.location.strand, end, distance, code, subject_start, False)
	
	if(potstart < potfinish):
		# Build the potential new feature for this location
		potfeat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(potstart, potfinish), strand = feat.location.strand)
		
		# Extract the sequence for this potential new feature
		potseq = SeqRecord.SeqRecord(potfeat.extract(seq_record.seq))
		
		# Count the stops in this sequence and return true if less than or equal to 1
		return([stopcount(potseq, table, 1, t) for t in [True, False]])
	else:
		return([100, 100])

def check_consistent_frame(stopdata):
	#stopdata = data
	return(len({d[1] for d in stopdata}) == 1)

def get_current_results(results, stopdata):
	locs = {d[0] for d in stopdata}
	return({i:l for i, l in results.items() if i in locs})

def find_and_filter_frame(currresults, feat, code, end, distance, subject_start, seq_record, table):
	#currresults, table = [results, translation_table]
	
	if(len(currresults) == 0):
		return(currresults, 1)
	
	# Generate list of lists containing [location, frame in subject, nstops inc final, nstops not inc final] 
	data = [[i, i%3 + 1] + get_stopcounts(i, l, feat, code, end, distance, subject_start, seq_record, table) for i, l in currresults.items()]
	
	# If any have one or zero stops, at the end, return these if OK
	checkdata = [d for d in data if d[2] <= 1 and d[3] == 0]
	if(len(checkdata) > 0):
		if(check_consistent_frame(checkdata) or len(checkdata) < 2):
			return(get_current_results(currresults, checkdata), 0)
		else:
			data = checkdata
	else:
		# Otherwise, i.e. there are some with internal stops 
		
		# Filter out any with greater than minimum number of nstops not inc final
		min_nsnf = min([d[2] for d in data])
		data = [d for d in data if d[2] == min_nsnf]
		if(check_consistent_frame(data) or len(data) < 2): return(get_current_results(currresults, data), 0)
		
		# Filter out any with greater than minimum number of nstops inc final
		min_nsif = min([d[3] for d in data])
		data = [d for d in data if d[3] == min_nsif]
		if(check_consistent_frame(data) or len(data) < 2): return(get_current_results(currresults, data), 0)
	
	## Filter out significantly uncommon frames
	#framecounts = Counter([d[1] for d in data])
	#maxfc = max([c for c in framecounts.values()])
	#data = [d for d in data if framecounts[d[1]] >= maxfc - 1]
	#if(check_consistent_frame(data) or len(data) < 2): return(get_current_results(currresults, data), 0)
	
	# Filter out frames where the nstops inc final do not match the stop status of the original feature (only if the original feature had no internal stops)
	endstop, instop = [stopcount(SeqRecord.SeqRecord(feat.extract(seq_record.seq)), table, 1, t) for t in [True, False]]
	if(instop == 0): data = [d for d in data if d[2] == endstop]
	if(check_consistent_frame(data) or len(data) < 2): return(get_current_results(currresults, data), 0)
	
	# Otherwise, return the results as they stand along with an error
	return(get_current_results(currresults, data), 1)

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

def correct_feature_by_query(feat, query_spec, seq_record, seqname, distance, featurename, translation_table, prioritise_longer_finishes):
	#feat, query_spec, distance, featurename, translation_table, prioritise_longer_finishes = [currfeat, stringspec, args.search_distance, name, args.translation_table, prioritise_long_stops]
	feat_start, feat_finish = feat.location.start, feat.location.end
	errstart = "Warning: sequence " + seqname + " has "
	
	distances_moved = {'start': 0, 'finish': 0}
	codon_start = None
	
	for end in ['start','finish']:
		#end = 'start'
		#end = 'finish'
		
		if(end not in query_spec.keys()):
			continue
		
		# Unpack search tuple
		code, query, out_rf, selector = query_spec[end] if len(query_spec[end]) == 4 else query_spec[end] + tuple("X")
		selector = out_rf if code == 'A' else selector
		
		errend = query  + " at the " + end + " of " + featurename + "\n"
		
		# Check if already ends with the searched sequence - REMOVED AS EXISTING SEQUENCE MAY BE SHORTER SUBSET OF DESIRED SEQUENCE e.g. TA TAA
		#if(end_already_correct(feat.extract(seq_record.seq), query, end, code, out_rf, args.translation_table)):
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
			target = [1,2,3,1,2][int(out_rf) + codon_start - 2] if(codon_start) else int(out_rf)
			
			if(end == "start"):
				results = {i:l for i, l in results.items() if (i - distance) % 3 + 1 == target }
				# Retain location if that location's rf (l+1)%3 is equal to the (target rf converted to subject rf)
			else:
				results = {i:l for i, l in results.items() if (abs(feat.location.end - feat_start) - distance + i - 1) % 3 + 1 == target}
			
			errmid = "no matches in the specified frame of " if len(results) == 0 else errmid
		
		truncated = False
		
		if(end == "start"):
			codon_start = 1
			
			# First check if truncated
			start_distance = len(seq_record) - feat.location.end if feat.location.strand == -1 else feat.location.start
			truncated = start_distance < distance
			
			# Retain only start locations that generate realistic amino acid sequences
			results, fail = find_and_filter_frame(results, feat, code, end, distance, subject_start, seq_record, translation_table)
			
			# Remove all results if completely impossible to determin
			if(fail): results = {}
			
			# If no feasible results at the start position, instead find the closest in-frame position 
			if(len(results) == 0 or truncated):
				errmid = "a truncated start or no ORF-producing matches (will set to closest ORF) of "
				
				# Set the current position - if normal, this is the current start position, if truncated, this is the end of the contig
				contig_start = list(find_all(str(subject_sequence), 'N'))[-1] + 1 if truncated else None
				current_position = contig_start if truncated else distance
				
				# Set the correction - if normal, this is 0, if truncated, this is 1, to ensure no results outside the contig
				correction = 1 if truncated else 0
				
				# Set up the three alternative results
				results = { current_position + v + correction : 1 for v in [-1, 0, 1] }
				
				# Find the result with suitable ORF
				results, fail = find_and_filter_frame(results, feat, code, end, distance, subject_start, seq_record, translation_table)
				if(fail): errmid = "no single detectable frame (no change will be made) in "
				
				# If truncated, set result to contig start but note codon position
				if(len(results) > 0 and truncated):
						codon_start = sorted(results.keys())[0] - contig_start + 1
						results = {contig_start : 1}
				
			
		else:
			
			if(len(results) > 0):
				# Prioritise longer matches 
				max_length = max(results.values())
				if(prioritise_longer_finishes): #Remove any matches shorter than the longest match
					results = {i:l for i, l in results.items() if l == max_length}
				elif(selector == "C"): # Remove matches shorter than the longest (or if longest >= 3 all matches) after the longest match
					first_max_i = min([i for i, l in results.items() if l == max_length])
					results = {i:l for i, l in results.items() if i <= first_max_i or (max_length < 3 and l == max_length)}
			
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
					location = locations[loc_dist.index(min(loc_dist))]
					if(loc_dist.count(min(loc_dist)) > 1): errmid = "multiple closest matches (taking first) of "
			
			
			feat_start, feat_finish = get_newends(location, results[location], feat_start, feat_finish, feat.location.strand, end, distance, code, subject_start, truncated)
			distances_moved[end] = location - distance
		else:
			errmid = "no succesful matches of " if errmid == "" else errmid
			codon_start = None
		
		if(errmid != ""):
			sys.stderr.write(errstart + errmid + errend)
	
	outfeat = copy.deepcopy(feat)
	if(feat_start < feat_finish):
		outfeat.location = SeqFeature.FeatureLocation(feat_start, feat_finish, feat.location.strand)
	else:
		distances_moved = {'start': 0, 'finish': 0}
	return(outfeat, codon_start, distances_moved)



def correct_truncated_features(feature, contiglength):
	#feature, contiglength = [feat, len(seq_record)]
	
	corrected_start, corrected_finish = [feature.location.start, feature.location.end]
	
	if(int(feature.location.start) == 0):
		corrected_start = SeqFeature.BeforePosition(corrected_start)
	
	if(int(feature.location.end) == contiglength):
		corrected_finish = SeqFeature.AfterPosition(corrected_finish)
	
	feature.location = SeqFeature.FeatureLocation(corrected_start, corrected_finish, feature.location.strand)
	

def end_already_correct(nuc_seq, query_seq, end, code, frame, translation_table):
	#nuc_seq = feat.extract(seq_record.seq)
	#query_seq = query
	#frame = out_rf
	
	
	if(code == 'N'):
		return(any(end == "start" and frame == "1" and nuc_seq.startswith(q) or # Starts with sequence, rf is 1
			   (end == "finish" and nuc_seq.endswith(q) and nuc_seq.rfind(q) % 3 + 1 == int(frame) )) for q in query_seq.split("/"))
	else:
		aa_seq = nuc_seq.translate(table = translation_table)
		return(any((end == "start" and aa_seq.startswith(q)) or
		           (end == "finish" and aa_seq.endswith(query_seq)) for q in query_seq.split("/")))

def check_context_features(contexts_in, seqname):
	#contexts_in = context_features
	# Check context_features all match in positions
	contexts_out = dict()
	
	for name, feats in contexts_in.items():
		#name, feats = list(context_features.items())[0]
		locations = [feat.location for feat in feats]
		
		feat_uniq = [feats[0]]
		for i in range(1, len(locations)):
			if(locations[i] != locations[i-1]):
				feat_uniq.append(feats[i])
		
		contexts_out[name] = feat_uniq
		
		if(len(feat_uniq) > 1):
			
			err = "Warning, positions of " + str(len(locations)) + " annotations for " + str(name) + " in " + str(seqname) + " do not match. If these are multiple distinct loci, this should be fine, but if these should cover the same locus, you may get incorrect results.\n"
			
			for i, feat in enumerate(feats):
				err += "\t(" + str(i+1) + ") " + feat.type +" is located at bases " + str(int(feat.location.start)+1) + " to " + str(int(feat.location.end)) + "\n"
		
			sys.stderr.write(err)
	
	return(contexts_out)
