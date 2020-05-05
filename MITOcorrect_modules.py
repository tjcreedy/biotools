#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 09:46:36 2020

@author: thomas
"""

import sys
import urllib.request
import re
import os
import csv
import subprocess

from collections import defaultdict #, Counter
from math import ceil, copysign
from Bio import Seq, SeqFeature, SeqRecord, SeqIO, AlignIO
from Bio.Align import AlignInfo

# =============================================================================
# 
# import copy
# from statistics import mode #, stdev
# import random, string
# import random, string
# 
# from Bio.Align import AlignInfo
# 
# from Bio import BiopythonWarning
# import warnings
# with warnings.catch_warnings():
#     warnings.simplefilter('ignore', BiopythonWarning)
# 
# =============================================================================
 
def loadnamevariants(source=None):
    variants = {}
    types = {}
    # Identify source
    
    if(source is None):
        url = 'https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt'
        source = urllib.request.urlopen(url)
    else:
        source = open(source, 'r')
    
    #Read source
     
    for line in source:
        line = line.decode('utf-8').strip()
        name = line.split(";")[0]
        annotype = line.split(":")[0].split(";")[1]
        types[name] = annotype
        listvars = line.split(":")[1].split(",")
        for v in listvars:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['',' GENE', ' '+annotype.upper()]:
                    variants[v+s] = name
     
    # Close handle
    source.close()
    # TODO check this works!
    return(variants, types)


def str_is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def parsespecs(path, alignpath, namevariants):
    #path, alignpath = [args.specifications, args.alignmentpaths]
    
    # Parse the master specifications file
    
    specs = dict()
    
    sys.stdout.write("Parsing specifications file %s\n" % (path))
    
    with open(path, 'r') as sh:
        #sh = open(path, 'r')
        
        # Extract the header row and split to list
        header = sh.readline().strip().split('\t')
        
        # TODO: check for duplicate headers
        
        ln = 1
        for line in sh:
            #line = sh.readline()
            
            ln += 1
            
            # Extract the values of the row and split to list
            items = line.strip().split('\t')
            
            # Get correct gene name
            name = None
            if items[0].upper() in namevariants:
                name = namevariants[items[0].upper()]
            
            if name is None:
                sys.exit(("Error: gene name %s in first column of line %s is ", 
                          "not recognised") % (items[0], str(ln)) )
            
            # Check start/stop column
            if items[1] not in ['start', 'stop']:
                sys.exit(("Error: end specification %s in second column of ",
                          "line %s is not recognised") % (items[1], str(ln)) )
            end = items[1]
            
            # Generate holding dict for this line's specs
            
            hold = dict()
            
            # Work through specifications
            for spec, value in zip(header[2:], items[2:]):
                #spec, value = [header[2], items[2]]
                
                # Do not record if value is empty
                if value == '':continue
                
                # If a list specification, generate list
                if '/' in value: 
                    value = value.split('/')
                
                # If a dict specification, generate dict
                if (type(value) is list and ',' in value[0] or 
                    type(value) is str and ',' in value):
                    
                    if (type(value) is list):
                        value = [v.split(',') for v in value]
                    else:
                        value = [value.split(',')]
                    
                    value = {k:v for k, v in value}
                
                # Add to specs dict
                
                hold[spec] = value
                
            
            # Do input checking
            for spec, value in hold.items():
                #spec = 'alignbody'
                #value = hold[spec]
                if spec == 'overlap':
                    # Overlap should be dict of recognised context names and 
                    # integer distances
                    if type(value) is not dict:
                        sys.exit(("Error: specification for %s on line %s is ",
                                 "not recognised") % (spec, str(ln)))
                    
                    for c, d in value.items():
                        #c, d = list(value.items())[0]
                        del hold[spec][c]
                        
                        if c in namevariants:
                            c = namevariants[c]
                        else: 
                            sys.exit(("Error: context name %s on line %s is ",
                                     "not recognised") % (c, str(ln)))
                        if str_is_int(d):
                            d = int(d)
                        else:
                            sys.exit(("Error: overlap %s for name %s on line ",
                                      "%s is not an integer") % 
                                     (d, c, str(ln)))
                        hold[spec][c] = d
                    
                elif spec in ['overlapmaxdistance', 'searchdistance']:
                    # MaxContextDistance, SearchDistance, 
                    # should all be > 0 integers
                    if not str_is_int(value):
                        sys.exit("Error: value %s for %s on line %s is not an integer" %
                                 ( value, spec, str(ln) ))
                    value = int(value)
                    if value < 1:
                        sys.exit("Error: value %s for %s on line %s is not greater than 0" %
                                 ( value, spec, str(ln) ))
                    
                    hold[spec] = value
                elif spec == 'alignbody':
                    # AlignBody should be a positive integer for start and a 
                    # negative integer for stop
                    if (str_is_int(value) and (
                            (end == 'start' and int(value) >= 0) 
                            or (end == 'stop' and int(value) <= 0))):
                        hold[spec] = int(value)
                    else:
                        sys.exit(("Error: value %s for %s on line %s is not",
                                  "zero or a %s integer") % 
                                 (value, spec, str(ln), 
                                  'positive' if end == 'start' else 'negative')
                                 )
                    
                elif spec == 'searchcode':
                    # Search code should be a single character A or N
                        if value not in 'AN' or len(value) != 1:
                            sys.exit("Error: value %s for %s on line %s is not A or N" %
                                     ( value, spec, str(ln) ))
                     
                elif spec == 'searchsequence':
                    # SearchSequence should be list of characters from one of
                    # two sets of possibilities for A or N respectively
                    if 'searchcode' not in hold:
                        sys.exit("Error: searchcode is required for searchsequence")
                    if type(value) is not list:
                        value = [value]
                        hold[spec] = value
                    allchar = list(''.join(value))
                    testchar = 'GPAVLIMCFYWHKRQNEDST*' 
                    if hold['searchcode'] == 'A': testchar = 'ATGC'
                    fails = [c for c in allchar if c not in testchar]
                    if len(fails) > 0:
                        sys.exit("Error: character(s) %s for %s on line %s are not valid" % 
                                 ( ', '.join(fails), spec, str(ln) ))
                        
                elif spec == 'searchreadframe':
                    # SearchReadFrame should be 1, 2 ,3 or *, and is irrelevant
                    # if using SearchCode A
                    if 'searchcode' not in hold:
                        sys.exit("Error: searchcode is required for searchsequence")
                    if value not in '123*' or len(value) > 1:
                        sys.exit("Error: value %s for %s on line %s is not valid" %
                                 ( value, spec, str(ln) ))
                    if hold['searchcode'] == 'A':
                        sys.stderr.write("Warning: searchreadframe is redundant for searchcode A on line %s")
                    if str_is_int(value):
                        hold[spec] = int(value)
                        
                elif spec == 'searchmultiselect':
                    # SearchMultiSelect should be one of five specific strings
                    if value not in ['F', 'FC', 'C', 'LC', 'L'] :
                        sys.exit("Error: value %s for %s on line %s is not valid" %
                                 ( value, spec, str(ln) ))
                        
                else:
                    sys.exit("Error: column %s is not recognised" %
                             (spec))
            
            # Initialise subdict if not already
            if(name not in specs): specs[name] = dict()
            
            # Input holding dict into main dict
            specs[name][end] = hold
            
    
    # Parse the alignment files if present
    
    if alignpath is not None:
        sys.stdout.write("Parsing alignment file %s\n" % (alignpath))
        
        with open(alignpath, 'r') as sh:
            #sh = open(alignpath, 'r')
            
            ln = 0
            for line in sh:
                #line = sh.readline()
                
                ln += 1
                
                # Extract the values of the row and split to list
                items = line.strip().split('\t')
                
                # Get correct gene name
                name = None
                if items[0].upper() in namevariants:
                    name = namevariants[items[0].upper()]
                
                if name is None:
                    sys.exit("Error: gene name %s on line %s is not recognised" %
                             ( items[0], str(ln) ))
                    
                elif name not in specs:
                    sys.exit("Error: gene name %s on line %s is not in %s" %
                             ( name, str(ln), path))
                
                specs[name]['apath'] = items[1]
    
    # TODO: come back to this?
# =============================================================================
#     # Generate a guide of what has full specifications, i.e. is to be done
#     
#     guide = dict()
#     
#     # List the required specifications for each category
#     speccat = {'overlap':   ['overlapmaxdistance',
#                              'overlap'],
#                'search':    ['searchdistance',
#                              'searchcode',
#                              'searchsequence',
#                              'searchmultiselect'],
#                'alignment': ['alignmentsearchdistance',
#                              'alignmentsegment',
#                              'apath']
#                }
#     
#     for target in specs:
#         guide[target] = dict()
#         # For each category, find the number of required specs
#         for cat, names in speccat.items():
#             n = len(names)
#             guideval = []
#             for e, v in specs[target]:
#                 # Count the number present for this target and end
#                 speccount = len([s for s in v if s in names])
#                 # If there are no specifications or if the number equals or 
#                 # exceeds the required number, record if this target/end is
#                 # to be processed, otherwise error out
#                 if speccount == 0 or speccount >= n:
#                     guideval.extend(speccount >= n)
#                 else:
#                     sys.exit("Error: %s specifications for %s %s are incomplete")
#             guide[target][cat] = any(guideval)
# =============================================================================
    return(specs)

def start_statsfile(outdir):
    
    # Open file
    stats = open(os.path.join(outdir, "filtering_results.tsv"), 'w')
    statwrite = csv.writer(stats, delimiter = '\t', quotechar = '', 
                               quoting = csv.QUOTE_NONE, escapechar = '')
    # Write header
    statwrite.writerow((['file', 'sequence_name', 'annotation_name',
                        'start_position', 'end_position', 'length', 
                        'start_match', 'end_match', 'reading_frame', 
                        'reading_frame_relative_to_original', 
                        'internal_stop_count']
                        + ["%s_%s_distance" % (t, e)
                           for e in ['start', 'end'] 
                           for t in ['overlap', 'consensus', 'body']]
                        + ['consensus_agreements', 'deletions', 
                           'insertions', 'final_score', 'selected']
                        ))
    
    # Return the handle and writer object
    return(stats, statwrite)



def get_file_details(filepath):
    sys.stdout.write("Starting work on genbank file %s\n" % (filepath))
    handle = SeqIO.parse(filepath, 'genbank')
    name = os.path.basename(filepath) 
    return(handle, name, 0)

def write_genbank_file(seqrecords, outdir, filename):
    outfile = os.path.join(outdir, filename)
    SeqIO.write(seqrecords, outfile, 'genbank')
    return([])

def get_features(seqrecord, namevariants):
    
    # Set up output containers
    features = defaultdict(list)
    unidentifiable_features = []
    unrecognised_names = set()
    other_features = []
    
    # Set up information variables
    nametags = ['gene', 'product', 'label', 'standard_name']
    othertypes = ['source', 'misc_feature', 'repeat_region', 'D-loop', 
                  'rep_origin','gap']
    
    for feat in seqrecord.features:
        #feat = seqrecord.features[1]
        # Remove any translations
        if('translation' in feat.qualifiers.keys()):
            del(feat.qualifiers['translation'])
        
        # Extract the tag that contains the feature name
        featname = None
        
        if(feat.type in othertypes):
            other_features.append(feat)
            continue
        elif(any(t in feat.qualifiers.keys() for t in nametags)):
            for t in nametags:
                if( t in feat.qualifiers.keys()):
                    featname = feat.qualifiers[t][0].upper()
                    break
        else:
            unidentifiable_features.append(feat)
            continue
        
        # Find the standard name
        if(featname in namevariants):
            features[namevariants[featname]].append(feat)
        else:
            unrecognised_names.add(featname)
    
    # Log
    unfeat = len(unidentifiable_features) > 0
    log = '1. Undentifiable features'
    if unfeat:
        log += '\n'
        for f in unidentifiable_features:
            log += "\t\t%s: %s-%sbp\n" % (f.type, str(f.location.start),
                                          str(f.location.end))
    else:
        log += ': NONE\n'
    return(features, unrecognised_names, unfeat,
           other_features + unidentifiable_features, log)

def clean_features(features, types):
    #features, types = feats, annotypes
    clean = dict()
    
    for name, feats in features.items():
        #name, feats = list(features.items())[-1]
        targets = [f for f in feats if f.type in ['CDS', 'tRNA', 'rRNA']]
        genes = [f for f in feats if f.type == 'gene']
        
        if len(targets) > 0:
            # TODO: remove duplicates if share same location details?
            
            clean[name] = set(targets)
        else:
            for gene in genes:
                gene.type = types[name]
            clean[name] = set(genes)
    
    return(clean)

def check_targets(clean, targets):
    #clean, targets, log =[cleanfeats, specs.keys(), log]
    absent, present, duplicate = [[], [], []]
    log = ''
    
    for target in targets:
        if target in clean:
            if len(clean[target]) == 1:
                present.append(target)
            else:
                duplicate.append(target)
        else:
            absent.append(target)
    
    for n, t, i in zip([2, 3], ['Absent', 'Duplicate'], [absent, duplicate]):
        log += "%s. %s target features" % (str(n), t)
        if(len(i) > 0):
            log += "\n\t%s\n" % (', '.join(i))
        else:
            log += ': NONE\n'
    
    return(present, log)

def gtl(lis):
    return(lis[0] > lis[1])

def ltl(lis):
    return(lis[0] < lis[1])

def overlap(initpos, strand, feats, specs, seqrecord):
    #strand, specs, feats = [feat.location.strand, specifications[target], cleanfeats]
    # Set up general variables
    absent = set()
    overdist = set()
    changes = [0, 0]
    
    # Work through the ends
    for i, end in enumerate(specs.keys()):
        # i, end = list(enumerate(specs.keys()))[1]
        # Extract the relevant specifications
        snames = ['overlapmaxdistance', 'overlap']
        maxdist, cspecs = [0, 0]
        if all([s in specs[end] for s in snames]):
            maxdist, cspecs = [specs[end][s] for s in snames]
        else:
            continue
        
        # Set up list for results
        correction = []
        
        # Generate a list of lists of all specified context features
        
        contexts = []
        for cname, dist in cspecs.items():
            #cname, dist = list(cspecs.items())[0]
            # Check if present, if not skip to next
            if cname not in feats:
                absent.add(cname)
                continue
            
            for cfeat in feats[cname]:
                contexts.append([cname, cfeat, dist])
        
        # Work through the specified context annotations
        for context in contexts:
            # context = contexts[1]
            cname, cfeat, specdist = context
            
            # Extract the feature locations
            cpos = [int(cfeat.location.start), int(cfeat.location.end)]
            
            # Find gap between context and target
            #gap = [min(cpos) - max(initpos), min(initpos) - max(cpos)]
            #gap = [g for g in gap if g > 0]
            
            # Set orientation (+ve, context follows target)
            orientation = 0
            # Set current distance between selected positions
            distance = 0
            # Set index of closest position
            closi = None
            
            # Find cross-wise distances
            crossdist = [max(cpos) - min(initpos), min(cpos) - max(initpos)]
            
            # Find the structure of the overlap
            if(abs(crossdist[0]) == abs(crossdist[1])):
                continue
            elif(abs(crossdist[1]) < abs(crossdist[0])):
                # Target is before context, overlap should be latter position 
                # of target and first position of context
                orientation = 1
                distance = crossdist[1]
                closi = initpos.index(max(initpos))
            else:
                # Target is after context, overlap should be first position 
                # of target and latter position of context
                orientation = -1
                distance = crossdist[0]
                closi = initpos.index(min(initpos))
            
            # Check gap and index of closest position are permissible
            
            if abs(distance) > maxdist + specdist:
                overdist.add(cname)
                continue
            
            if closi != i:
                continue
            
            # Calculate the exact new position
            change = distance + orientation * specdist
            newpos =  initpos[i] + change
            
            # If the position is suitable, store it
            if (0 <= newpos <= len(seqrecord) and          # inside the contig
                # On +ve strand, newstart is less than (before) current stop 
                # or on -ve stand, current stop is before newstart
                ((i == 0 and ltl([newpos, initpos[1]][::strand])) or
                 # On +ve strand, newstop is greater than (after) current start
                 # or on -ve strand, current start is after newstop
                 (i == 1 and gtl([newpos, initpos[0]][::strand])))):
                 
                correction.append([cname, newpos, change * strand])
        
        if len(correction) > 0:
            cval = [abs(cv[2]) for cv in correction]
            mini = cval.index(min(cval))
            changes[i] = correction[mini]
        
    
    logstring = []
    output = [0, 0]
    for i, e in enumerate(['start', 'stop']):
        c = changes[i]
        if type(c) is list:
            logstring.append("%sbp from %s according to %s" %
                      (str(c[2]), e, c[0]))
            output[i] = c[1]
        else:
            logstring.append("0bp from %s" % (e))
            output[i] = initpos[i]
    if len(absent) > 0:
        logstring.append("%s not present" % (', '.join(absent)) )
    
    log = '\t\tOverlap positions: %s\n' % (', '.join(logstring))
    
    return(output, log)

def relative_rf(target, reference, referencerf, strand):
    return(( strand * (target - reference) % 3 + referencerf) % 3)

def repeat_series(series, length, si=0):
    denom = len(series)
    out = [series[i % denom] for i in range(length+si)]
    return(out[si:length+si])

def get_regions(initpos, adjpos, rf, specs, strand, seqrecord, table):
    #adjpos, specs, strand, table = [contextpos, specifications[target], feat.location.strand, args.translationtable]
    
    #Set up output
    output = dict()
    
    for i, end in enumerate(['start', 'stop']):
        #i, end = list(enumerate(specs.keys()))[1]
        
        # Extract the distance
        dist = specs[end]['searchdistance']
        
        # Convert distance for amino acids, ensuring complete codons
        dist = ceil(dist/3) * 9 if specs[end]['searchcode'] == 'A' else dist
        
        # Delimit the region
        regpos = [adjpos[i] - dist, adjpos[i] + dist]
        
        # Truncate if exceeds the contig
        regpos[0] = 0 if regpos[0] < 0 else regpos[0]
        regpos[1] = len(seqrecord) if regpos[1] > len(seqrecord) else regpos[1]
        
        # Flip if strand is reverse
        regpos = regpos[::strand]
        
        
        # Truncate if exceeds the other end of the annotation
        # If end is start and, on positive strand, the existing end is before 
        # (less than) the search end, or on strand is negative, the search end
        # is before the existing end, set the search end to the existing end
        if i == 0 and ltl([adjpos[1], regpos[1]][::strand]):
            regpos[1] == adjpos[1]
        # And vice versa
        if i == 1 and gtl([adjpos[0], regpos[0]][::strand]):
            regpos[0] == adjpos [0]
        
        # Find length
        reglen = max(regpos) - min(regpos)
        
        # Get the frame of the first base in the region:
        regrf = 0
        regrf = relative_rf(regpos[0], initpos[0], rf, strand)
        
        # Generate the feature
        subject_feat = SeqFeature.SeqFeature(
                SeqFeature.FeatureLocation(min(regpos), max(regpos)), 
                strand=strand)
        sequence = subject_feat.extract(seqrecord.seq)
        if specs[end]['searchcode'] == 'A':
            sequence = sequence.translate(table=table)
        
        rfcodes = [0, 1, 2]
        # SEQquence, POSition of each base in contig, Position Codon Position
        # of each base relative to annotation reading frame
        output[end] = {'seq': str(sequence),
                       'pos': list(range(min(regpos), max(regpos)+1)
                                   )[::strand][:-1],
                       'pcp': repeat_series(rfcodes, reglen, 
                                            rfcodes.index(regrf))}
    
    return(output, '')


def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def get_search_results(regions, adjpos, specs, seqrecord, strand, table, ff):
    # regions, adjpos, specs, strand, table, ff = [searchregions, contextpos, specifications[target], feat.location.strand, args.translationtable, args.framefree]
    
    results = dict()
    truncated = ''
    
    # Work through ends and specs
    for end in ['start', 'stop']:
        #end = 'stop'
        sspecs = specs[end]
        
        # Set default codon_start for hits
        cs = 0
        
        # Find locations of query relative to query frame and positions
        localresults = dict()
        for q in sspecs['searchsequence']:
            #q = sspecs['searchsequence'][0]
            # Work through locations of hits
            for i in find_all(regions[end]['seq'], q):
                #i = list(find_all(regions[end]['seq'], q))[0]
                # If the current hit location exists and is longer than the
                # current hit, do not change, else add current hit length
                if (i not in localresults or 
                    (i in localresults and len(q) > localresults[i][1])):
                    localresults[i] = [cs, len(q)]
        
        # If region appears to be truncated at the start, add further results
        if end == 'start' and regions[end]['pos'][0] in [0, len(seqrecord)]:
            localresults.update( {0: [i, 1] for i in [0,1,2]} )
            truncated = 'start'
        
        # If region appears to be truncated at the stop, remove short matches
        # and add a match at the end
        if (not ff and end == 'stop' 
            and regions[end]['pos'][-1] in [0, len(seqrecord)]
            and sspecs['searchcode'] is 'N'):
            
            localresults = {i: l for i, l in localresults.items if l[1] >= 3}
            localresults[regions[end]['pos'][-1]] = 1
            truncated = 'stop' if not truncated else 'startstop'
        
        # TODO: add in default locations around existing initial positions?
        
        # Convert locations to contig positions and annotation frame
        outresults = dict()
        for i, l in localresults.items():
            # Dict where key is position of start of annotation in contig
            outresults[regions[end]['pos'][i]] = {
                    # li = position in search region, to allow for correct
                    # ordering of positions if strand is reverse
                    'li': i,
                    # len = The length of the match
                    'len': l[1],
                    # The codon_start of the match, 0-indexed
                    'cs': l[0],
                    # The reading frame of the match relative to the original
                    # annotation, modified by the codon_start so that start and
                    # stop matches in the correct frame are correctly paired
                    'arf': (regions[end]['pcp'][i] + l[0]) % 3}
        
        
        # If a stop region, remove all matches after the first >= 3 match if 
        # the largest match is >=3, on a frame-by-frame basis, unless the 
        # search sequences are not inframe (i.e. freeframe is on)
        if not ff and end == 'stop' and not truncated:
            for prf in [0, 1, 2]:
                #prf = 1
                res_prf = ({k: v for k, v in outresults.items() 
                           if v['arf'] == prf })
                if len(res_prf) > 0:
                    maxmatch = max([v['len'] for v in res_prf.values()])
                    if maxmatch >= 3:
                        firsti = min([v['li'] for v in res_prf.values() if
                                      v['len'] >= 3])
                        for k, v in list(outresults.items()):
                            if v['li'] > firsti and v['arf'] == prf:
                                del outresults[k]
        
        results[end] = outresults
    
    # Generate the sequence records for all combinations
    options = []
    for startpos, startdet in results['start'].items():
        for stoppos, stopdet in results['stop'].items():
            #startpos, startdet = list(results['start'].items())[-1]
            #stoppos, stopdet = list(results['stop'].items())[1]
            
            # Correct the stop position by the length
            stoppos = stoppos + stopdet['len'] * strand
                # This shouldn't be theoretically possible, but just in case...
            if stoppos < 0: stoppos = 0
            if stoppos > len(seqrecord): stoppos = len(seqrecord)
            
            # Set up a list of start and stop positions
            pos = [startpos, stoppos]
            
            # Skip the combination if:
                # the start position is not less than the stop position (or vice
                # versa on reverse strand)
            if (gtl(pos[::strand])
                # the total length of the resulting sequence is less than the
                # sum of the lengths of the matches
                or max(pos) - min(pos) < startdet['len'] + stopdet['len']
                # the results are not in the same reading frame unless search
                # sequences are not expected to be in frame (i.e. freeframe on)
                or (not ff and startdet['arf'] != stopdet['arf'])):
                continue
            
            # Set the positions to be truncated if necessary
            for i, e in enumerate(['start', 'stop']):
                # For each truncated end, apply the appropriate position type
                if e in truncated:
                    # Set as 
                    # BeforePosition if start and forward or stop and reverse
                    # AfterPosition if stop and forward or start and reverse
                    pos[i] = (SeqFeature.BeforePosition(pos[i]) 
                              if strand == [1, -1][i]
                              else SeqFeature.AfterPosition(pos[i]))
            
            
            # Calculate the distance from the overlap position
            # TODO: Check correct in -ve strand - should be!
            overlapdist = [p -a for a, p in zip(adjpos, pos)]
            
            # Make the feature
            pos = pos[::strand]
            newfeat = SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(pos[0], pos[1], strand=strand))
            newfeat.qualifiers['codon_position'] = startdet['cs'] + 1
            
            # Extract the nucleotide sequence
            ntseq = newfeat.extract(seqrecord.seq)
            
            # Check the length and end
            trail = len(ntseq) % 3
            
            # Pad the sequence for amino acid extraction
            transseq = ntseq
            
            if trail > 0:
                trailbases = str(ntseq[-trail:])
                pad = 'A' if trailbases == ['T','TA'][trail-1] else 'N'
                transseq += Seq.Seq(pad * (3 - trail))
            
            # Extract the amino acids
            aaseq = transseq.translate(table = table)
            
            # Count stops
            aastring = str(aaseq)
            internal = re.sub("\**$", '', aastring)
            intstops = internal.count('*')
            termstops = re.findall("\**$", aastring)[0]
            termstops = 0 if len(termstops) == 0 else len(termstops[0])
            
            # Gather data and append
            options.append({'arf': startdet['arf'], 'feat': newfeat,
                            'lens': [d['len'] for d in [startdet, stopdet]],
                            'inst': intstops, 'tmst': termstops,
                            'nt': ntseq, 'aa': aaseq,
                            'adjd': overlapdist})
            
            # TODO some logging here?
    
    return(options, '')

def filter_searchresults(results):
    #results = searchresults
    
    # Filter out any with more than the minimum number of internal stops
    intstops = [r['inst'] for r in results]
    minintstop = min(intstops)
    minintstop = 0 if minintstop > 3 else minintstop
    results = [r for r in results if r['inst'] == minintstop]
    
    # Filter out those with the lowest overlap score. Always retain a minimum
    # number of results, plus a proportion of the remainder, removing the rest
    overlapscore = [sum([abs(a) for a in r['adjd']]) for r in results]
    retentionthresh = 6
    retentionprop = 0.51
    if (len(overlapscore) > retentionthresh + round(retentionprop)):
        nretain = round(retentionprop * (len(overlapscore) - retentionthresh))
        maxval = sorted(overlapscore)[retentionthresh + nretain - 1]
        results = [r for r, s in zip(results, overlapscore) if s <= maxval]
    
    # TODO: more filtering? Perhaps rf based?
    
    # TODO some logging here?
    return(results, '')


def ungapped_distance(seq, value, ref):
    # seq, value, ref =  [ntalign, seqss[e], conss[e]]
    locs = [value, ref]
    between = seq[min(locs):max(locs)]
    dist = len(between) - between.count('-')
    return(copysign(1, value-ref) * dist)

def cumsum(lis):
    total = 0
    for i in lis:
        total += i
        yield total

def gapped_distance(seq, value, dist):
    # seq, value, dist = [consensus, conss['stop'], stdist['stop']]
    seqv = [0 if s == '-' else 1 for s in seq]
    out = None
    if dist < 0:
        seqv = [(-1 * l) + 1 for l in cumsum(seqv[:value:][::-1])]
        out = value - seqv.index(dist)
    else:
        seqv = [c - 1 for c in cumsum(seqv[value:])]
        out = value + seqv.index(dist)
    if out < 0: out = None
    if out > len(seq): out = None
    
    return(out)

def align_and_analyse(results, args, specs, target, seqname, temp):
    # results, specs, aligntype = [filterresults, specifications[target], args.alignmenttype]
    
    # Align and generate alignment scores
    for result in results:
        #result = results[0]
        
        # Set up the filenames for input and output
        
        name = "%s_%s_%s-%s" % (seqname, target, 
                                int(result['feat'].location.start) + 1,
                                int(result['feat'].location.end))
        suffixes = ['_result.fasta', '_align.fasta']
        files = [os.path.join(temp, name + s) for s in suffixes]
        
        # Generate a temporary fasta
        outrecord = SeqRecord.SeqRecord(seq=result[args.alignmenttype], 
                                        id=name, description='')
        SeqIO.write(outrecord, files[0], 'fasta')
        
        # Build and run the mafft command:
            # Generate the commands
        arguments = ("mafft --quiet --6merpair --addfragments".split(' ') 
                     + [files[0], specs['apath']])
            # Open the output file for writing
        with open(files[1], 'w') as out:
            # Run the process
            p = subprocess.Popen(arguments,
                             stdout = out,
                             universal_newlines = True)
            # Wait to finish
            p.communicate()
        
        # Read the alignment
        alignment = AlignIO.read(files[1], 'fasta')
        
        # Find consensus
        alignment_summary = AlignInfo.SummaryInfo(alignment)
        consensus = str(alignment_summary.gap_consensus())
        
        # Find target sequence
        ntalign = str(alignment[-1].seq)
        
        # Delete the files
        os.remove(files[0])
        if not args.keepalignments: os.remove(files[1])
        
        # Find start and stop positions 
        seqss = dict()
        conss = dict()
        
        for e, r in zip(['start', 'stop'], ["^-*", "-*$"]):
            #e, r = list(zip(['start', 'stop'], ["^-*", "-*$"]))[0]
            s, c = [len(re.search(r, seq).group()) 
                    for seq in [ntalign, consensus]]
            if e == 'stop':
                s, c = [len(consensus) - v for v in [s, c]]
            seqss[e], conss[e] = [s, c]
        
        # Generate standard body distances
        stdist = {e: specs[e]['alignbody'] for e in ['start', 'stop']}
        #stdist = {'start': 0, 'stop': -4} # COX1
        bodyss = {e: gapped_distance(consensus, conss[e], v) for 
                  e, v in stdist.items()}
        
        # Determine distances
        cond = [0, 0]
        bodd = [0, 0]
        
        for i, e in enumerate(['start', 'stop']):
            # i, e = list(enumerate(['start', 'stop']))[1]
            if (e == 'start' and seqss[e] > conss[e] or 
                e == 'stop' and seqss[e] < conss[e]):
                cond[i] = ungapped_distance(consensus, seqss[e], conss[e])
            else:
                cond[i] = ungapped_distance(ntalign, seqss[e], conss[e])
            
            if (e == 'start' and seqss[e] > bodyss[e] or 
                e == 'stop' and seqss[e] < bodyss[e]):
                bodd[i] = (ungapped_distance(consensus, seqss[e], bodyss[e])
                             - stdist[e])
            else:
                bodd[i] = (ungapped_distance(ntalign, seqss[e], bodyss[e])
                             + stdist[e])
        
        diserr = cond, bodd
        if args.alignmenttype == 'AA':
            diserr = [[i * 3 for i in j] for j in diserr]
        
        
        # Determine agreement, deletion or insertion for each of the positions
        # of the ntseq after trimming terminal gaps
        cid = []
        for n, c in zip(ntalign[seqss['start']:seqss['stop']], 
                        consensus[seqss['start']:seqss['stop']]):
            if n == c or (n != '-' and c != '-'):
                cid.append('-')
            elif n == '-' and c != '-':
                cid.append('d')
            elif n != '-' and c == '-':
                cid.append('i')
        
        ccts = [cid.count(c) for c in '-di']
        
        
        # Generate a score for each end of the annotation
        endscore = [None, None]
        for i in [0, 1]:
            
            # TODO - step through this manually because it's giving strange
            # values...
            
            # Generate the weighted average of the alignment distances
            alignscore = (abs(ccts[i]) * args.alignmentweight
                          + abs(bodd[i]) * (1 - args.alignmentweight)) / 2
            
            # Generate the weighted average of the overlap dist and alignscore
            endscore[i] = (abs(result['adjd'][i]) * args.overlapweight
                           + alignscore * (1 - args.overlapweight)) / 2
        
         # Compile
        result.update({'cond': cond, 'bodd': bodd, 'ccts': ccts, 
                       'slen': len(result['nt']), 'score': sum(endscore) / 2})
    
    # Find the minimum score and mark the selection in the results
    minscore = min([r['score'] for r in results])
    for r in results:
        r['select'] = r['score'] == minscore
    
    return(results, '')

def write_detailed_results(results, gbname, seqname, target, sw=None):
    # results = alignresults
    # Start output list
    
    statl = []
    for result in results:
        # Extract the feature
        feat = result['feat']
        
        # Construct a line to write to the stats file
        stats = ([gbname, seqname, target,
                  str(feat.location.start + 1),
                  str(feat.location.end),
                  result['slen'],
                  str(result['nt'][:result['lens'][0]]),
                  str(result['nt'][-result['lens'][1]:]),
                  feat.qualifiers['codon_position'],
                  str(result['arf'] + 1),
                  result['inst']]
                  + [result[k][i] for i in [0, 1] 
                   for k in ['adjd', 'cond', 'bodd']]
                  + result['ccts'] 
                  + [result['score'], result['select']])
        if sw:
            sw.write(stats)
        else:
            statl.append(stats)
    
    return(statl)

def generate_output_target(results, target, args):
    
    resultfeats = list()
    
    for result in results:
        if args.potentialfeatures:
            
            #result = results[0]
            
            # Extract the feature
            feat = result['feat']
            
            # Generate name for this result
            name = "%s_%s-%s" % (target, int(feat.location.start) + 1,
                                 int(feat.location.end))
            
            # Give the feature the name and a type
            feat.qualifiers['gene'] = name
            feat.qualifiers['note'] = 'potential_MITOcorrect_annotation'
            feat.type = 'potential'
            
            # Add the feature to the outputs
            resultfeats.append(feat)
            
        elif result['select']:
            
            feat = result['feat']
            
            # Add feat details
            feat.qualifiers['gene'] = target
            feat.type = 'CDS'
            
            # TODO: add second feat for gene
            
            # Add the feature to the outputs
            resultfeats.append(feat)
            
    
    return(resultfeats)


def initialise(args):
    
    # Read in the namevariants file
    namevariants, annotypes = loadnamevariants()
    
    # Allow user to pass an additional namevariants file
    if args.namevariants:
        morevariants, moretypes = loadnamevariants(args.namevariants)
        namevariants.update(morevariants)
        annotypes.update(moretypes)
    
    # Read in and parse specifications table to a dict 
    specs = parsespecs(args.specifications, args.alignmentpaths, 
                                    namevariants)
    
    # Make output directory and temporary directory, open logfile
    temp = os.path.join(args.outputdirectory, 'intermediate_alignments')
    if not os.path.exists(temp):
        os.makedirs(temp)
    
    log = (open(os.path.join(args.logfile), 'w') 
        if args.logfile else open(os.devnull, 'w'))
    
    # Start an output file for writing statistics if args.detailedresults
    if args.detailedresults:
        stath, statw = start_statsfile(args.outputdirectory)
    
    return(namevariants, annotypes, specs, temp, log, stath, statw)

def prepare_seqrecord(seqn, seqrecord, gbname, namevariants, annotypes,
                      specifications):
    
    issues = defaultdict(set)
    
    log = '##### Sequence %s from file %s #####\n' % (seqrecord.name, gbname)
    
    sys.stdout.write("Correcting sequence %s: %s\n" % (seqn, seqrecord.name))
    
    # Extract target feature types from the record, other features, and 
    # variables for any issues encountered. Record any issues to the log
    
    features, unnames, unfeat, ofeats, flog = get_features(seqrecord, 
                                                          namevariants)
    
    if len(unnames) > 0: issues['unrecnames'].add(unnames)
    if unfeat: issues['hasunidfeats'].add(seqrecord.name)
    
    # Clean the features to reject any annotations of the non-target type
    
    cleanfeats = clean_features(features, annotypes)
    
    # Establish which of the features specified by the user are present
    present, clog = check_targets(cleanfeats, specifications.keys())
    log +=  flog + clog + '4. Target features'
    if len(present) == 0:
         log += 'NONE\n'
    
    # Append the non-target feats to the ofeats
    ofeats.extend([c for n, cs in cleanfeats.items() 
                   for c in list(cs) 
                   if n not in present])
    
    return(present, cleanfeats, ofeats, issues, log)

def correct_feature(cleanfeats, specifications, gbname, seqrecord, args, 
                    temp, target):
    # specifications, target = [specs, present[6]]
    
    # TODO: something to ensure original annotation is always part of the
    # results list
    
    # Extract the CDS
    feat = list(cleanfeats[target])[0]
    
    # Set up log and stats
    log = "\n\t%s\n" % (target) 
    statl = []
    
    # Determine the start and stop positions of the current annotation,
    # slicing by strand (1 or -1) makes order correct
    initpos = [int(feat.location.start),
               int(feat.location.end)][::feat.location.strand]
    
    # Determine the reading frame of the current annotation if present
    # Reading frame is 0 indexed internally
    rf = 0
    if 'codon_position' in feat.qualifiers:
        cs = feat.qualifiers['codon_position']
        cs = cs[0] if type(rf) is list else cs
        rf = [0, 2, 1][cs-1]
    
    # Get context adjustments
    contextpos, clog = overlap(initpos, feat.location.strand, cleanfeats,
                               specifications[target], seqrecord)
    
    # Determine and extract the search regions
    searchregions, rlog = get_regions(initpos, contextpos, rf,
                                      specifications[target], 
                                      feat.location.strand,
                                      seqrecord, args.translationtable)
    
    # Find possible start and stop positions
    searchresults, slog = get_search_results(searchregions,
                                             contextpos,
                                             specifications[target],
                                             seqrecord,
                                             feat.location.strand,
                                             args.translationtable,
                                             args.framefree)
    
    # Make an empty list, the default output if no result
    result = []
    flog, alog = ['', '']
    if len(searchresults) > 0:
        
        # Filter results
        filterresults, flog = filter_searchresults(searchresults)
        
        # Align the filtered results and generate alignment stats
        if len(filterresults) > 0:
            alignresults, alog = align_and_analyse(filterresults, args,
                                                   specifications[target],
                                                   target, seqrecord.name, 
                                                   temp)
        else:
            alignresults, alog = [filterresults, '']
        
        if args.detailedresults:
            # Output a list of features to write to the contig
            # Write a table of filter results
            statl = write_detailed_results(alignresults, gbname,
                                                   seqrecord.name, target)
        
        # Output final result(s)
        # TODO: selection of single result if equal scores
        result = generate_output_target(alignresults, target, args)
        
        
        
        if args.potentialfeatures or len(result) == 0:
            result.append(feat)
    
    
    if len(result) == 0:
         # TODO: generate list of two features, target type and gene based on input feat
        result = [feat]
    
    
    # Extend outfeats with these
    return(result, [log + clog + rlog + slog + flog + alog], statl)