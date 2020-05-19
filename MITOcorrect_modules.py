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
import functools
import time
import datetime
from collections import defaultdict
from math import ceil, copysign
from Bio import Seq, SeqFeature, SeqRecord, SeqIO, AlignIO
from Bio.Align import AlignInfo

def elapsed_time(start):
    end = time.perf_counter()
    elapsed = round(end - start, 2)
    return(" | %ss | " % str(elapsed))

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
    
    sh = open(path, 'r')
    
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
        
    
    sh.close()
    # Parse the alignment files if present
    
    if alignpath is not None:
        sys.stdout.write("Parsing alignment file %s\n" % (alignpath))
        
        sh = open(alignpath, 'r')
        
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
        sh.close()
        
    return(specs)

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
    log = 'Undentifiable features: '
    if unfeat:
        loguf = []
        for f in unidentifiable_features:
            loguf.append("%s %s-%sbp" % (f.type, str(f.location.start),
                                           str(f.location.end)))
        log += ', '.join(loguf) + "\n"
    else:
        log += 'NONE\n'
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
    log = []
    
    for target in targets:
        if target in clean:
            if len(clean[target]) == 1:
                present.append(target)
            else:
                duplicate.append(target)
        else:
            absent.append(target)
    
    for t, i in zip(['Absent', 'Duplicate'], [absent, duplicate]):
        r = ', '.join(i) if len(i) > 0 else 'NONE'
        log.append("%s target features: %s\n" % (t, r))
    
    return(present, log)

def gtl(lis):
    return(lis[0] > lis[1])

def ltl(lis):
    return(lis[0] < lis[1])

def overlap(initpos, strand, feats, specs, seqrecord):
    #strand, specs, feats = [feat.location.strand, specifications[target], cleanfeats]
    # Set up general variables
    absent = defaultdict(set)
    overdist = defaultdict(set)
    changes = [None, None]
    
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
                absent[end].add(cname)
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
                overdist[end].add(cname)
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
    outpos = [0, 0]
    outchange = [None, None]
    for i, e in enumerate(['start', 'stop']):
        c = changes[i]
        mlog = ''
        if c:
            mlog += ("%sbp from %s according to %s" % (str(c[2]), e, c[0]))
            outpos[i] = c[1]
            outchange[i] = c[2]
        else:
            mlog += ("unchanged at %s" % (e))
            outpos[i] = initpos[i]
        elog = []
        if len(absent[e]) > 0:
            elog.append("%s not present" % (', '.join(absent[e])))
        if len(overdist[e]) > 0:
            elog.append("%s over distance" % (', '.join(absent[e])))
        if len(elog) > 0:
            mlog += (" (%s)" % ('; '.join(elog)))
        logstring.append(mlog)
    
    log = 'Overlap positions: %s\n' % (', '.join(logstring))
    
    return(outpos, outchange, log)

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
    log = dict()
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
        log[end] = "%s %sbp-%sbp" % (end, str(min(regpos)), str(max(regpos)))
        rfcodes = [0, 1, 2]
        # SEQquence, POSition of each base in contig, Position Codon Position
        # of each base relative to annotation reading frame
        output[end] = {'seq': str(sequence),
                       'pos': list(range(min(regpos), max(regpos)+1)
                                   )[::strand][:-1],
                       'pcp': repeat_series(rfcodes, reglen, 
                                            rfcodes.index(regrf))}
    log = ', '.join([log[e] for e in ['start', 'stop']])
    log = "Search regions: %s\n" % (log)
    return(output, log)


def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def get_search_results(regions, adjustment, specs, seqrecord, strand, table, ff):
    # regions, adjustment, specs, strand, table, ff = [searchregions, (contextpos, change), specifications[target], feat.location.strand, args.translationtable, args.framefree]
    # Adjustment is tuple containing adjpos and change
    
    results = dict()
    trunc = ''
    
    # Work through ends and specs
    for end in ['start', 'stop']:
        #end = 'start'
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
        
        # Convert localresults to list
        localresults = [[i]+v for i, v in localresults.items()]
        
        # Define the posible positions for a search region to start or stop on
        # in order to determine the feature to be truncated. Note these include
        # 1 position within the contig to deal with indexes not matching
        truncpos = [0, 1, len(seqrecord) - 1, len(seqrecord)]
        
        # If region appears to be truncated at the start, add further results
        if end == 'start' and regions[end]['pos'][0] in truncpos:
            localresults += [[0, i, 1] for i in [0, 1, 2]]
            trunc = 'start is'
        
        # If region appears to be trunc at the stop, remove short matches
        # and add a match at the end. Truncation is defined by the last position
        if (not ff and end == 'stop' 
            and regions[end]['pos'][-1] in truncpos
            and sspecs['searchcode'] is 'N'):
            
            localresults = [[i, c, l] for i, c, l in localresults if l >= 3]
            endpos = len(regions[end]['seq'])-1
            localresults += [[endpos, i, 1] for i in [0, 1, 2]]
            trunc = 'stop is' if not trunc else 'both start and stop are'
        
        # TODO: add in default locations around existing initial positions?
        
        # Convert locations to contig positions and annotation frame
        outresults = []
        for i, c, l in localresults:
            # Dict where key is position of start of annotation in contig
            outresults.append([regions[end]['pos'][i], {
                    # li = position in search region, to allow for correct
                    # ordering of positions if strand is reverse
                    'li': i,
                    # len = The length of the match
                    'len': l,
                    # The codon_start of the match, 0-indexed
                    'cs': c,
                    # The reading frame of the match relative to the original
                    # annotation, modified by the codon_start so that start and
                    # stop matches in the correct frame are correctly paired
                    'arf': (regions[end]['pcp'][i] + c) % 3}])
        
        
        # If a stop region, remove all matches after the first >= 3 match if 
        # the largest match is >=3, on a frame-by-frame basis, unless the 
        # search sequences are not inframe (i.e. freeframe is on)
        if not ff and end == 'stop' and not 'stop' in trunc:
            fixresults = []
            for prf in [0, 1, 2]:
                #prf = 0
                res_prf = [[k, v] for k, v in outresults if v['arf'] == prf ]
                if len(res_prf) > 0:
                    maxmatch = max([v['len'] for k, v in res_prf])
                    if maxmatch >= 3:
                        firsti = min([v['li'] for k, v in res_prf if 
                                      v['len'] >= 3])
                        for k, v in outresults:
                            if v['li'] <= firsti and v['arf'] == prf:
                                fixresults.append([k,v])
                    else:
                        fixresults.extend(res_prf)
            outresults = fixresults
        
        results[end] = outresults
    
    # Log so far
    log = ["%s hits at %s" % (str(len(results[e])), e) 
           for e in ['start', 'stop']]
    log = "Query search: found %s" % (' and '.join(log))
    
    # Generate the sequence records for all combinations
    options = []
    for startpos, startdet in results['start']:
        for stoppos, stopdet in results['stop']:
            #startpos, startdet = results['start'][1]
            #stoppos, stopdet = results['stop'][-1]
            
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
                if e in trunc:
                    # Set as 
                    # BeforePosition if start and forward or stop and reverse
                    # AfterPosition if stop and forward or start and reverse
                    pos[i] = (SeqFeature.BeforePosition(pos[i]) 
                              if strand == [1, -1][i]
                              else SeqFeature.AfterPosition(pos[i]))
            
            # Calculate the distance from the overlap position as long as the 
            # overlap position is a true overlap adjustment - if there has been
            # no movement, return 0 in all cases
            adjpos, change = adjustment
            overlapdist = [p - a if c is not None else 0 
                           for a, p, c in zip(adjpos, pos, change)]
            
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
                            'adjd': overlapdist, 'reject': None})
            
    
    log += ", generating %s possible result regions" % (str(len(options)))
    log += '; ' + trunc + ' truncated\n' if trunc else '\n'
    return(options, log)

def filter_searchresults(results):
    #results = searchresults
    resultsremain = len(results)
    # Filter out any with more than the minimum number of internal stops
    intstops = [r['inst'] for r in results]
    minintstop = min(intstops)
    maxis = 0 if minintstop > 3 else minintstop
    intstoprejects = 0
    for r in results:
        if r['inst'] > maxis:
            r['reject'] = "intstop > " + str(minintstop)
            intstoprejects += 1
            resultsremain -=1
    log = ("Filtering: %s results removed for > %s internal stops, " % 
           (str(intstoprejects), str(minintstop)))
    
    # Filter out those with the worst overlap score. Always retain a minimum
    # number of results, plus a proportion of the remainder, up to a maximum
#    overlapscore = [sum([abs(a) for a in r['adjd']]) for r in results]
#    retentionthresh = 5
#    retentionprop = 0.51
#    retentionmax = 20
#    overlapscorerejects = 0
#    if resultsremain > retentionthresh + round(retentionprop):
#        # Calculate the number to retain
#        noverthresh = resultsremain - retentionthresh
#        noverthreshretain = round(retentionprop * noverthresh)
#        nretain = retentionthresh + noverthreshretain
#        nretain = retentionmax if nretain > retentionmax else nretain
#        # Find the maximum overlap score to retain
#        maxval = sorted(overlapscore)[nretain - 1]
#        for r, s in zip(results, overlapscore):
#            if s > maxval and not r['reject']:
#                r['reject'] = "overlapscore >" + str(maxval)
#                overlapscorerejects += 1
#                resultsremain -= 1
#        log += ("%s results removed for > %s overlap score, " % (
#                                                      str(overlapscorerejects),
#                                                      str(maxval)))
#    # TODO: more filtering? Perhaps rf based?
#    
    log += "%s results remain \n" % (str(resultsremain))
    return(results, log)


def ungapped_distance(seq, value, ref):
    # seq, value, ref =  [targetalign, seqss[e], conss[e]]
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
    # results, specs, seqname, aligntype = [filterresults, specifications[target], seqrecord.name, args.alignmenttype]
    
    # Align and generate alignment scores
    for result in results:
        #result = results[0]
        # Fill with blank stats if already rejected
        if result['reject']:
            result.update({'cond': ['',''], 'bodd': ['',''],
                           'ccts': ['','',''], 'slen': len(result['nt']),
                           'score': ''})
            continue
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
            # TODO: add --anysymbol if using amino acid alignments
        arguments = "mafft --quiet --6merpair --addfragments"
        arguments = arguments.split(' ') + [files[0]]
        if args.alignmenttype == 'aa': arguments += ['--anysymbol'] 
        arguments += [specs['apath']]
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
        targetalign = str(alignment[-1].seq)
        # Delete the files
        os.remove(files[0])
        if not args.keepalignments: os.remove(files[1])
        # Find start and stop positions 
        seqss = dict()
        conss = dict()
        # Search for start and stop positions of the sequence and consensus
        for e, r in zip(['start', 'stop'], ["^-*", "-*$"]):
            #e, r = list(zip(['start', 'stop'], ["^-*", "-*$"]))[0]
            s, c = [len(re.search(r, seq).group()) 
                    for seq in [targetalign, consensus]]
            if e == 'stop':
                s, c = [len(consensus) - v - 1 for v in [s, c]]
            seqss[e], conss[e] = [s, c]
        # Generate standard body distances
        stdist = {e: specs[e]['alignbody'] for e in ['start', 'stop']}
        if args.alignmenttype == 'aa':
            stdist = {e: round(v/3) for e, v in stdist.items()}
        bodyss = {e: gapped_distance(consensus, conss[e], v) for 
                  e, v in stdist.items()}
        # Determine distances
        cond = [0, 0]
        bodd = [0, 0]
        #Find the distances for the sequence
        for i, e in enumerate(['start', 'stop']):
            # i, e = list(enumerate(['start', 'stop']))[1]
            if (e == 'start' and seqss[e] > conss[e] or 
                e == 'stop' and seqss[e] < conss[e]):
                cond[i] = ungapped_distance(consensus, seqss[e], conss[e])
            else:
                cond[i] = ungapped_distance(targetalign, seqss[e], conss[e])
            if (e == 'start' and seqss[e] > bodyss[e] or 
                e == 'stop' and seqss[e] < bodyss[e]):
                bodd[i] = (ungapped_distance(consensus, seqss[e], bodyss[e])
                             - stdist[e])
            else:
                bodd[i] = (ungapped_distance(targetalign, seqss[e], bodyss[e])
                             + stdist[e])
        # cond, bodd = two lists, each containing start and stop error for the 
        # consensus and the body respectively
        if args.alignmenttype == 'aa':
            cond, bodd = [[i * 3 for i in j] for j in [cond, bodd]]
        # Generate a score for each end of the annotation
        endscore = [None, None]
        for i in [0, 1]:
            
            # Generate the weighted average of the alignment distances
            alignscore = (abs(cond[i]) * args.alignmentweight
                          + abs(bodd[i]) * (1 - args.alignmentweight))
            
            # Generate the weighted average of the overlap dist and alignscore
            endscore[i] = (abs(result['adjd'][i]) * args.overlapweight
                           + alignscore * (1 - args.overlapweight))
            
        # Determine agreement, deletion or insertion for each of the positions
        # of the ntseq after trimming terminal gaps
        cid = []
        for n, c in zip(targetalign[seqss['start']:seqss['stop']], 
                        consensus[seqss['start']:seqss['stop']]):
            if n == c or (n != '-' and c != '-'):
                cid.append('-')
            elif n == '-' and c != '-':
                cid.append('d')
            elif n != '-' and c == '-':
                cid.append('i')
        # Generate list holding counts of agreements, deletions and insertions
        ccts = [cid.count(c) for c in '-di']
        # Compile
        result.update({'cond': cond, 'bodd': bodd, 'ccts': ccts, 
                       'slen': len(result['nt']),
                       'score': round(sum(endscore) / 2, 2)})
    
    # Find the minimum score and mark the selection in the results
    scores = [r['score'] for r in results if r['score'] is not '']
    minscore = min(scores) if len(scores) > 0 else ''
    selected = []
    for r in results:
        if r['score'] is not '':
            if r['score'] > minscore:
                out = 'score > ' + str(round(minscore,2))
                r['reject'] = out
            else:
                selected.append("result region %sbp-%sbp" %
                                    (str(int(r['feat'].location.start) + 1),
                                     str(int(r['feat'].location.end))))
    log = "Alignment: "
    success = False
    if len(scores) < 1:
        log += "no results remain to align and analyse\n"
    elif len(selected) < 1:
        log += "selected no regions from %s results\n" % (str(len(scores)))
    else:
        success = True
        log += "selected %s of %s" % (', '.join(selected), str(len(scores)))
        log += " results with minimum combined overlap and alignment score of "
        log += "%s\n" % (str(minscore))
    
    return(results, success, log)

def write_detailed_results(results, gbname, seqname, target):
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
                  + [result['score'], result['reject']])
        statl.append(stats)
    
    return(statl)

def generate_output_target(results, target, args):
    #results = alignresults
    
    resultfeats = []
    
    for result in results:
        
        if args.potentialfeatures:
            
            #result = results[0]
            
            # Extract the feature
            feat = result['feat']
            
            # Generate outcome
            outcome = result['reject'] if result['reject'] else "selected"
            # Generate name for this result
            name = "%s_%s-%s: %s" % (target, int(feat.location.start) + 1,
                                 int(feat.location.end), outcome)
            
            # Give the feature the name and a type
            feat.qualifiers['gene'] = name
            feat.qualifiers['note'] = 'potential_MITOcorrect_annotation'
            feat.type = 'potential'
            
            # Add the feature to the outputs
            resultfeats.append(feat)
            
        elif not result['reject']:
            
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
    
    sys.stdout.write("Completed initialisation, starting processing\n")
    
    return((namevariants, annotypes, specs, temp))

def prepare_seqrecord(seqrecord, gbname, namevariants, annotypes, 
                      specifications, pid, logq):
    start = time.perf_counter()
    issues = defaultdict(set)
    
    log = 'PID%s file %s sequence %s' % (pid, gbname, seqrecord.name)
    
    # Extract target feature types from the record, other features, and 
    # variables for any issues encountered. Record any issues to the log
    
    features, unnames, unfeat, ofeats, flog = get_features(seqrecord, 
                                                          namevariants)
    logq.put(log + flog)
    if len(unnames) > 0: issues['unrecnames'].add(unnames)
    if unfeat: issues['hasunidfeats'].add(seqrecord.name)
    
    # Clean the features to reject any annotations of the non-target type
    
    cleanfeats = clean_features(features, annotypes)
    
    # Establish which of the features specified by the user are present
    present, clog = check_targets(cleanfeats, specifications.keys())
    for c in clog:
        logq.put(log + c)
    
    # Generate target log
    tlog =  'Target features: '
    if len(present) == 0:
        tlog += 'NONE\n'
    else:
        tlog += "%s present to process\n" % (str(len(present)))
    
    # Append the non-target feats to the ofeats
    ofeats.extend([c for n, cs in cleanfeats.items() 
                   for c in list(cs) 
                   if n not in present])
    
    logq.put(log + elapsed_time(start) + tlog)
    return(present, cleanfeats, ofeats, issues)

def correct_feature(cleanfeats, specifications, gbname, seqrecord, args, 
                    temp, pid, logq, target):
    
    # specifications, target = [specs, present[0]]
    # specifications, target = [specs, 'CYTB']
    
    # TODO: something to ensure original annotation is always part of the
    # results list
    
    # Extract the CDS
    feat = list(cleanfeats[target])[0]
    
    # Set up log and stats
    statl = []
    log = "PID%s %s %s" % (pid, seqrecord.name, target)
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
    start = time.perf_counter()
    contextpos, change, clog = overlap(initpos, feat.location.strand,
                                       cleanfeats, specifications[target],
                                       seqrecord)
    logq.put(log + elapsed_time(start) + clog)
    # Determine and extract the search regions
    start = time.perf_counter()
    searchregions, rlog = get_regions(initpos, contextpos, rf,
                                      specifications[target], 
                                      feat.location.strand,
                                      seqrecord, args.translationtable)
    logq.put(log + elapsed_time(start) + rlog)
    # Find possible start and stop positions
    start = time.perf_counter()
    searchresults, slog = get_search_results(searchregions,
                                             (contextpos, change),
                                             specifications[target],
                                             seqrecord,
                                             feat.location.strand,
                                             args.translationtable,
                                             args.framefree)
    logq.put(log + elapsed_time(start) + slog)
    
    # Make an empty list, the default output if no result
    result = []
    if len(searchresults) > 0:
        # Filter results
        start = time.perf_counter()
        filterresults, flog = filter_searchresults(searchresults)
        logq.put(log + elapsed_time(start) + flog)
        # Align and generate alignment stats
        start = time.perf_counter()
        alignresults, alog = align_and_analyse(filterresults, args,
                                               specifications[target],
                                               target, seqrecord.name, temp)
        logq.put(log + elapsed_time(start) + alog)
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
    else:
         # TODO: generate list of two features, target type and gene based on input feat
        result = [feat]
    
    # Extend outfeats with these
    return(result, statl)

def get_seqrecords(filepaths, onefile):
    '''Generator outputting the file name and seqrecord for each entry in each
       input file'''
    # filepaths, onefile = [args.genbank, args.onefile]
    fpiter = iter(filepaths)
    count = 0
    currfile = next(fpiter, None)
    while currfile:
        # If expecting to output to one file overall, keep the count iterating
        # across input files; otherwise, reset the count for each input file
        if not onefile: count = 0
        # Check to see if the current file is the last one
        nextfile = next(fpiter, None)
        # Set up to iterate on the current input file
        sriter = SeqIO.parse(currfile, 'genbank')
        # Set the name for this file
        outname = onefile if onefile else currfile
        # Get the first record from this file
        currrecord = next(sriter, None)
        # Work through the records
        while currrecord:
            count +=1
            # Check to see if the current record is the last one by attempting
            # to get the next one
            nextrecord = next(sriter, None)
            # If we have reached the final seqrecord for the current file, and 
            # either we are outputting file by file or we are outputting in one
            # file and we've finished the last file, set the value of  count to
            # the total number of seqrecords to output for the relevant output 
            # file, otherwise 0 denoting we still don't know the total number
            filecount = 0 
            if not nextrecord and (not onefile or (onefile and not nextfile)):
                filecount = count
            # Yield a tuple containing the file name to output to, the record
            # and the count for this filename
            yield(currfile, outname, currrecord, filecount)
            # Move the next record to the current record for the next iteration
            currrecord = nextrecord
            
        # Move the next file to the current file for the next iteration
        currfile = nextfile

def write_stats(outdir, statq):
    # Open file handles
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
                           'insertions', 'final_score', 'outcome']
                        ))
    stats.flush()
    # Wait on items from queue and write them as received
    while 1:
        queueitem = statq.get()
        if queueitem == 'kill': break
        for l in queueitem:
            statwrite.writerow(l)
        stats.flush()
    stats.close()

def write_log(outdir, logfile, logq):
    # Set up the logging output handle
    logh = None
    if logfile:
        logh = open(os.path.join(outdir, logfile), 'w')
    # Start a constant process that waits to receive data in the form of a file
    # specifying the source of the seqrecord, and the new seqrecord
    while 1:
        logline = logq.get()
        if logline == 'kill': break
        if logfile:
            logh.write(logline)
            logh.flush()
    if logfile: logh.close()

def write_genbanks(outdir, filepaths, onefile, seqq):
    '''Sets up the output file writing handles then waits to receive data
       requesting printing of log info and seqrecord to the specified file'''
    #outdir, filepaths, onefile, queueitem = [args.outputdirectory, args.genbank, args.onefile, (outname, seqrecord, filetotal)]
    # Set up the gbhands dict to work whether all sequences are to be output
    # to the same file or different files.
    gbh = dict()
    if onefile:
        outh = open(os.path.join(outdir, onefile), 'w')
        gbh[onefile] = {'h': outh,
                        'c': 0,
                        't': 0}
    else:
        for file in filepaths:
            outpath = os.path.join(outdir, file)
            # For multiple files, for each file set up dict with handle, 
            # counter and total number of seqrecords to expect.
            gbh[file] = {'h': open(outpath, 'w'),
                         'c': 0,
                         't': 0}
    # Start a constant process that waits to receive data in the form of a file
    # specifying the source of the seqrecord, and the new seqrecord
    while 1:
        queueitem = seqq.get()
        if queueitem == 'kill': break
        file, seqrecord, filetotal = queueitem
        # Filetotal will be either 0 or the total number of seqrecords for
        # this file. Each time a seqrecord is received, increment the count
        # The total will only be the correct value once the one seqrecord
        # carrying the total arrives (the last one read), and the handle
        # will only close when the counter equals the total
        gbh[file]['h'].write(seqrecord.format("genbank"))
        gbh[file]['h'].flush()
        gbh[file]['c'] += 1
        gbh[file]['t'] += filetotal
        if gbh[file]['c'] == gbh[file]['t']:
            gbh[file]['h'].close()

def print_terminal(filenames, prinq):
    #filenames = args.genbank
    # Get total
    tot = 0
    for file in filenames:
        with open(file, 'r') as fh:
            for line in fh:
                if re.search('^LOCUS', line):
                    tot += 1
    # Set up to print
    start = time.perf_counter()
    done = 0
    remain = "unknown time"
    current = ''
    run = True
    # Print
    while run:
        sys.stdout.write("\r%s\r" % (' ' * 70))
        line = ("\rCorrected %s of %s total records, %s remaining%s" % 
                (done, tot, remain, current))
        sys.stdout.write(line)
        sys.stdout.flush()
        queueitem = prinq.get()
        if queueitem == 'kill': break
        done += 1
        now = time.perf_counter()
        elapsed = now-start
        remain = round((elapsed/done) * (tot - done))
        remain = "approx " + str(datetime.timedelta(seconds=remain))
        #current = ", currently processing " + queueitem
    now = time.perf_counter()
    elapsed = round(now-start)
    elapsed = str(datetime.timedelta(seconds=elapsed))
    elapsedper = str(datetime.timedelta(seconds=elapsed/done))
    line = "\nFinished in %s, %s per record\n" % (elapsed, elapsedper)
    sys.stdout.write(line)
    sys.stdout.flush()
    run = False



#def print_terminal(filenames, prinq):
#    #filenames = args.genbank
#    # Count the number of expected seqrecords
#    ctr = dict()
#    for file in filenames:
#        ctr[file] = [0, 0, 0, 0, 0] # Done, total, proportion, pbarw, %
#        with open(file, 'r') as fh:
#            for line in fh:
#                ctr[file][1] += 1 if re.search('^LOCUS', line) else 0
#    # Calculate the total
#    tot = sum([v[1] for k, v in ctr.items()])
#    
#    
##    # Set up the counter
##    stdscr = curses.initscr()
##    curses.noecho()
##    curses.cbreak()
#    
#    start = time.perf_counter()
#    done = 0
##    barw = 15
#    while 1:
#        # Get the incoming completed file part and check
#        file = prinq.get()
#        if not file:
#            break
#        done += 1
#        # Calculate timings
#        now = time.perf_counter()
#        elapsed = round(now-start)
#        #elapstr = str(datetime.timedelta(seconds=elapsed))
#        remain = (elapsed/done) * (tot - done)
#        remainstr = str(datetime.timedelta(secords=remain))
##        # Iterate the revelant file counter
##        ctr[file][0] += 1
##        ctr[file][2] = (ctr[file][0]/ctr[file][1])
##        ctr[file][3] = round(ctr[file][2] * barw)
##        ctr[file][4] = round(ctr[file][2] * 100)
##        # Build the screen
##        line = "Total processed: {0}. Elapsed time: {1}. \
##                Approx time remaining: {3}".format(done, elapstr, remainstr)
##        stdscr.addstr(0, 0, line)
##        for i, f in enumerate(ctr.keys()):
##            l = "{0} : [{4:" + str(4 + barw) + "}] {1}/{2} ({3}%)".format(
##                    f, ctr[f][0], ctr[f][1], ctr[f][4], "#" * ctr[f][3])
##            stdscr.addstr(i + 1, 0, l)
##        stdscr.refresh()
#        sys.stdout.write("Done %s of %s, approx %s remaining\r" % (done, tot,
#                                                                   remainstr))
#        sys.stdout.flush()
#    curses.echo()
#    curses.nocbreak()
#    curses.endwin()


def start_writers(pool, manager, args):
    
    seqq = manager.Queue()
    statq = manager.Queue() if args.detailedresults else None
    logq = manager.Queue()
    prinq = manager.Queue()
    
    seqwatch = pool.apply_async(functools.partial(write_genbanks,
                                                  args.outputdirectory,
                                                  args.genbank,
                                                  args.onefile),
                                (seqq,))
    
    if args.detailedresults:
        statwatch = pool.apply_async(functools.partial(write_stats, 
                                                       args.outputdirectory),
                                     (statq,))
    
    logwatch = pool.apply_async(functools.partial(write_log,
                                                  args.outputdirectory,
                                                  args.logfile),
                                    (logq,))
    
    printwatch = pool.apply_async(functools.partial(print_terminal, 
                                                    args.genbank),
                                  (prinq,))
    
    return(seqq, statq, logq, prinq, 
           (seqwatch, statwatch, logwatch, printwatch))

def process_seqrecord(args, utilityvars, seqq, statq, logq, prinq, indata):
    #
    # indata = next(seqrecordgen)
    gbname, outname, seqrecord, filetotal = indata
    namevars, annotypes, specs, temp = utilityvars
    
    pid = os.getpid()
    # Extract the necessary items from the seqrecord and clean
    present, cleanfeats, ofeats, issues = prepare_seqrecord(seqrecord, gbname,
                                                            namevars, 
                                                            annotypes, specs, 
                                                            pid, logq)
    
    # Process the present cleanfeatures
    if len(present) > 0:
        outfeats = []  
        for target in present:
            outfeat, statl = correct_feature(cleanfeats, specs, gbname,
                                             seqrecord, args, temp, pid, 
                                             logq, target)
            outfeats.extend(outfeat)
            if args.detailedresults: statq.put(statl)
        
        # Replace all features with the new ones and add on the others
        if len(outfeats) > 0:
            seqrecord.features = outfeats + ofeats
    else:
        seqrecord.features = seqrecord.features
        # TODO log that no target features found and that the sequence is
        # being output as-is
    
    seqq.put((outname, seqrecord, filetotal))
    prinq.put(gbname)