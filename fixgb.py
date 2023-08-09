#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 09:52:39 2020

@author: thomas
"""

# Imports

import sys
import re
import os
import urllib
import argparse
import textwrap as _textwrap

import Bio
from Bio import SeqIO
from copy import copy
from collections import defaultdict


#TODO:
# Suppress biopython warnings by default - e.g.
# /usr/lib/python3/dist-packages/Bio/GenBank/__init__.py:1395: BiopythonParserWarning: Expected sequence length 5289, found 5125 (urn.local...13-b53ordm).
#   warnings.warn(
# Add argument to just read and write GB using biopython if no other biopython arguments are called - to fix lengths like the above
# Add code to replace ACCESSION, VERSION, DEFINITION if . or urn.local* or anything else unwanted (grep "ACCESSION" on all gb to check)
# Fix incorrect origin wrapping e.g.
# BIOD01157.gb
# /usr/lib/python3/dist-packages/Bio/GenBank/__init__.py:361: BiopythonParserWarning: Attempting to fix invalid location '29865..15146' as it looks like incorrect origin wrapping. Please fix input file, this could have unintended behavior.
#   warnings.warn(
# Traceback (most recent call last):
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/__init__.py", line 1126, in location
#     cur_feature.location = SeqFeature.FeatureLocation(
#   File "/usr/lib/python3/dist-packages/Bio/SeqFeature.py", line 811, in __init__
#     raise ValueError(
# ValueError: End location (15146) must be greater than or equal to start location (29864)
# During handling of the above exception, another exception occurred:
# Traceback (most recent call last):
#   File "/home/thomc/programming/bioinformatics/biotools//fixgb.py", line 354, in <module>
#     for i, seqrecord in enumerate(gbin):
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/Scanner.py", line 516, in parse_records
#     record = self.parse(handle, do_features)
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/Scanner.py", line 499, in parse
#     if self.feed(handle, consumer, do_features):
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/Scanner.py", line 470, in feed
#     self._feed_feature_table(consumer, self.parse_features(skip=False))
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/Scanner.py", line 420, in _feed_feature_table
#     consumer.location(location_string)
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/__init__.py", line 1131, in location
#     cur_feature.location = _loc(
#   File "/usr/lib/python3/dist-packages/Bio/GenBank/__init__.py", line 369, in _loc
#     f1 = SeqFeature.FeatureLocation(s_pos, expected_seq_length, strand)
#   File "/usr/lib/python3/dist-packages/Bio/SeqFeature.py", line 811, in __init__
#     raise ValueError(
# ValueError: End location (15997) must be greater than or equal to start location (29864)
# Is there a way to change the origin???
# Add option and code to fix IDs with _reversed or _modified


# Global variables

divisions = 'BCT,PRI,ROD,MAM,VRT,INV,PLN,VRL,PHG,RNA,SYN,UNA,EST,STS,GSS,HTG,UNK'.split(',')

# Class definitons


class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Monkey Patches


Bio.SeqFeature.FeatureLocation.__hash__ = lambda self : hash(f"{self.start.__str__()}"
                                                             f"{self.end.__str__()}"
                                                             f"{self.strand}")
Bio.SeqFeature.SimpleLocation = Bio.SeqFeature.FeatureLocation

# Function definitions


def loadnamevariants():
    conversion = {}
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    fullparse = {}
    alltypes = set()
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])

        fullvariants = []
        for v in [name] + variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' '+annotype.upper()]:
                    fullvariants.append(v+s)
                    conversion[v+s] = name

        alltypes.add(annotype)
        fullparse[name] = {'type': annotype, 'variants': fullvariants}
    return conversion, alltypes, fullparse


def checkn(lis, name, ln, alt = None):
    n = len(lis)
    if n == 0: return alt
    if n == 1: return lis[0]
    if n > 1: sys.exit("Line %s has multiple matches for %s" % (str(ln), name))

def fixhead(inpath, outpath, args):
    # Open file(s) and set up handler(s)
    ingb = open(inpath, 'r') if type(inpath) is str else inpath
    outgb = open(outpath, 'w') if type(outpath) is str else outpath

    ln = 0

    # Work through file
    opentag = False
    tagstart = "                     /"
    for line in ingb:
        # line = ingb.readline()
        line = line.rstrip()
        ln += 1
        # Check if header line
        if re.match("^LOCUS", line):

            # Length
            length = re.findall(r'\d+ +(?:bp(?= ))', line)
            length = checkn(length, "sequence length", ln, '  bp')

            # Type
            seqtype = [s for s in ['DNA', 'RNA', 'tRNA', 'mRNA', 'uRNA', 'snRNA', 'cDNA'] if
                       ' ' + s + ' ' in line]
            seqtype = '{:<7}'.format(checkn(seqtype, "sequence type", ln, ''))

            # Strand
            strand = [s for s in ['ss-', 'ds-', 'ms-'] if ' ' + s + ' ' in line]
            strand = checkn(strand, "sequence strand", ln, '   ')

            # Date
            date = re.findall(r'(?<= )\d{2}-[A-Z]{3}-\d{4}', line)
            date = checkn(date, "date", ln, '01-JAN-1979')

            # Division
            division = [s for s in divisions if ' ' + s + ' ' in line]
            division = checkn(division, "division", ln, args.defaultdivision)

            # Topology
            topology = [s for s in ['linear', 'circular'] if ' ' + s + ' ' in line]
            topology = '{:<8}'.format(checkn(topology, "topology", ln, ''))

            # Name
            remaining = line
            for i in ['LOCUS', length, seqtype, strand, date, division, topology]:
                if i is not None and not i.isspace():
                    remaining = remaining.replace(i.strip(), '')

            name = remaining.split()
            notname = [n for n in name if n != line.split()[1]]
            if len(notname) > 0:
                sys.exit(f"Line {ln} has unidentified value(s): {', '.join(notname)}")

            name = checkn(name, "sequence name", ln, None)
            if name is None: sys.exit(f"Line {ln} has no identifiable sequence name")

            # Format together
            # Do name and length
            length = length.split()[0]
            nllens = [len(name), len(length)]
            if sum(nllens) > 27:
                sys.exit(f"Line {ln} sequence name {name} is too long, must be {27 - nllens[1]} or "
                         f"fewer characters to allow length")
            namelength = f"{name}{' ' * (28 - sum(nllens))}{length}"

            # Output
            line = f"LOCUS{' ' * 7}{namelength} bp {strand}{seqtype} {topology} {division} {date}"

        # Fix open tags by removing line with open tag and starting next line with a /
        elif line == tagstart:
            opentag = True
            continue
        elif opentag:
            line = tagstart + line.lstrip()
            opentag = False

        outgb.write(f"{line}\n")

    if type(inpath) is str:
        ingb.close()

    if type(outpath) is str:
        outgb.close()

def sortextractfeats(record, extracttypes):
    focalfeats = {}
    otherfeats = []
    for feat in record.features:
        if feat.type in extracttypes:
            if feat.location in focalfeats:
                focalfeats[feat.location].append(feat)
            else:
                focalfeats[feat.location] = [feat]
        else:
            otherfeats.append(feat)
    return focalfeats, otherfeats


def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


def set_feat_name(feat, name):
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys() :
                feat.qualifiers[t][0] = name
    return feat


def overlap(feats):
    a, b = [[int(e) for e in (f.location.start, f.location.end)] for f in feats]
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]


def remove_duplicates(feats):
    featdict = defaultdict(list)

    for feat in feats:
        id = f"{get_feat_name(feat)} {feat.type}"
        featdict[id].append(feat)

    outfeats = []
    for id, featlis in featdict.items():
        featcheck = copy(featlis)
        i = 0
        while i < len(featcheck):
            overlapj = [j for j, f in enumerate(featcheck) if overlap([f, featcheck[i]])]
            if len(overlapj) > 1:
                lens = [len(featcheck[j]) for j in overlapj]
                maxj = [j for j, l in zip(overlapj, lens) if l == max(lens)][0]
                dropj = [j for j in overlapj if j != maxj]
                featcheck = [f for j, f in enumerate(featcheck) if j not in dropj]
                nregress = sum(1 for j in dropj if j <= i)
                i -= nregress
            i += 1
        outfeats.extend(featcheck)
    return outfeats


def add_anticodon(feats, nameconvert):
    # feats = trnas
    outfeats = []
    noanticodons = 0
    # Iterate
    for i, featlis in enumerate(feats.values()):
        # featlis = list(feats.values())[1]
        # Check if at least one is a tRNA, dicard all if not
        types = set(f.type for f in featlis)
        if 'tRNA' not in types:
            outfeats.extend(featlis)
            continue
        # Get the names of the features
        names = [get_feat_name(feat) for feat in featlis]
        trnanames = [get_feat_name(feat) for feat in featlis if feat.type == 'tRNA']
        #sys.stderr.write(f"{i} {names}\n")
        # Check if the names are already recognisable, discard any that are, discard any that don't
        # have the same name as any trnas
        featcont = []
        for feat, name in zip(featlis, names):
            if name in nameconvert or name not in trnanames:
                outfeats.append(feat)
            else:
                featcont.append(feat)
        if len(featcont) == 0:
            continue
        # Search for anticodon tags
        anticodons = {}
        for feat in featcont:
            # feat = featcont[1]
            if 'anticodon' in feat.qualifiers:
                name = get_feat_name(feat)
                # Check for tags
                acstring = feat.qualifiers['anticodon']
                if len(acstring) > 1:
                    sys.stderr.write(f"Warning: {feat.type} annotation {name} has multiple "
                                     f"anticodon tags, this will be skipped\n")
                    continue
                # Search for information within tag
                acstring = acstring[0]
                parser = {#'position': r"pos:([0-9]+)\.\.([0-9]+)",
                          'aminoacid': r"aa:([A-Za-z]+)",
                          'anticodon': r"seq:([A-Za-z]{3})"}
                acdata = {}
                for part, regex in parser.items():
                    m = re.search(regex, acstring)
                    if not m:
                        sys.stderr.write(f"Warning: {part} cannot be found in {feat.type} "
                                         f"annotation {name} anticodon tag, this will be skipped\n")
                        continue
                    acdata[part] = m.groups(1) if len(m.groups(1)) > 1 else m.groups(1)[0]
                if len(acdata) < 2:
                    continue
                # Construct putative new names
                acname = f"TRN{acdata['aminoacid'].upper()[0]}-{acdata['anticodon'].upper()}"
                if name in anticodons:
                    anticodons[name].add(acname)
                else:
                    anticodons[name] = {acname}
        if len(anticodons) == 0:
            noanticodons += 1
            outfeats.extend(featcont)
            continue
        # Set up renamer
        rename = {}
        for fname, acname in anticodons.items():
            if len(acname) > 1:
                sys.stderr.write(f"Warning: anticodon tags do not match for the two or more "
                                 f"annotations named {fname} at the same locus")
                outfeats.extend(featcont)
                continue
            rename[fname] = list(acname)[0]
        # Rename all
        outfeats.extend([set_feat_name(feat, rename[get_feat_name(feat)]) for feat in featcont])
    # Output list of features
    return outfeats, noanticodons


def add_anticodon(feats, nameconvert):
    # feats = trnas
    outfeats = []
    noanticodons = 0
    # Iterate
    for i, featlis in enumerate(feats.values()):
        # featlis = list(feats.values())[1]
        # Check if at least one is a tRNA, dicard all if not
        types = set(f.type for f in featlis)
        if 'tRNA' not in types:
            outfeats.extend(featlis)
            continue
        # Get the names of the features
        names = [get_feat_name(feat) for feat in featlis]
        trnanames = [get_feat_name(feat) for feat in featlis if feat.type == 'tRNA']
        #sys.stderr.write(f"{i} {names}\n")
        # Check if the names are already recognisable, discard any that are, discard any that don't
        # have the same name as any trnas
        featcont = []
        for feat, name in zip(featlis, names):
            if name in nameconvert or name not in trnanames:
                outfeats.append(feat)
            else:
                featcont.append(feat)
        if len(featcont) == 0:
            continue
        # Search for anticodon tags
        anticodons = {}
        for feat in featcont:
            # feat = featcont[1]
            if 'anticodon' in feat.qualifiers:
                name = get_feat_name(feat)
                # Check for tags
                acstring = feat.qualifiers['anticodon']
                if len(acstring) > 1:
                    sys.stderr.write(f"Warning: {feat.type} annotation {name} has multiple "
                                     f"anticodon tags, this will be skipped\n")
                    continue
                # Search for information within tag
                acstring = acstring[0]
                parser = {#'position': r"pos:([0-9]+)\.\.([0-9]+)",
                          'aminoacid': r"aa:([A-Za-z]+)",
                          'anticodon': r"seq:([A-Za-z]{3})"}
                acdata = {}
                for part, regex in parser.items():
                    m = re.search(regex, acstring)
                    if not m:
                        sys.stderr.write(f"Warning: {part} cannot be found in {feat.type} "
                                         f"annotation {name} anticodon tag, this will be skipped\n")
                        continue
                    acdata[part] = m.groups(1) if len(m.groups(1)) > 1 else m.groups(1)[0]
                if len(acdata) < 2:
                    continue
                # Construct putative new names
                acname = f"TRN{acdata['aminoacid'].upper()[0]}-{acdata['anticodon'].upper()}"
                if name in anticodons:
                    anticodons[name].add(acname)
                else:
                    anticodons[name] = {acname}
        if len(anticodons) == 0:
            noanticodons += 1
            outfeats.extend(featcont)
            continue
        # Set up renamer
        rename = {}
        for fname, acname in anticodons.items():
            if len(acname) > 1:
                sys.stderr.write(f"Warning: anticodon tags do not match for the two or more "
                                 f"annotations named {fname} at the same locus")
                outfeats.extend(featcont)
                continue
            rename[fname] = list(acname)[0]
        # Rename all
        outfeats.extend([set_feat_name(feat, rename[get_feat_name(feat)]) for feat in featcont])
    # Output list of features
    return outfeats, noanticodons

def getcliargs(arglist=None):

    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
    Standalone tool for fixing various issues in genbank files, supplied on STDIN
    |n
    To fix issues so the file can be opened with BioPython, run with --fixhead. The header will be
    checked for correct spacing, and missing division codes will be filled in with UNK by default, 
    unless a three-letter code is supplied to -d/--defaultdivision. Any incomplete annotation tags
    will be fixed as far as possible.
    |n
    On mitogenomes, fill missing annotations with --fillpairs. This will check through 
    the annotations to make sure that all 'CDS', 'tRNA' or 'rRNA' annotations have an identical 
    'gene' annotation, and vice-versa. All other types of annotation will be ignored. 
    |n
    Where multiple overlapping annotations share the same name and type, retain only the longest 
    annotation using --removeduplicates
    |n
    Standardise gene names with --standardname. This will check each gene name and replace it with
    the standard name used in https://github.com/tjcreedy/constants.
    |n
    Remove any geneious editing history or annotation derivation annotations using --removegeneious.
    |n
    By default, any source annotations will be reset to start at position 1 and end at the last 
    position in the sequence. Use --dontresetsource to stop this behaviour.
    """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("--fixhead", action='store_true',
                        help="fix issues with headers so they work in BioPython")
    parser.add_argument("-d", "--defaultdivision", default="UNK", metavar="DIV",
                        help="a three-letter code to be used as the default devision where the "
                             "division code is missing")
    parser.add_argument("--fixlengths", action='store_true',
                        help="repair length values in LOCUS line")
    parser.add_argument("--addanticodons", action='store_true',
                        help="try to add the anticodon to the name of tRNAs without it")
    parser.add_argument("--removeduplicates", action='store_true',
                        help="remove shorter duplicate annotations")
    # parser.add_argument("-g", "--deletegeneious", action='store_true',
    #                     help="try to remove automatic Geneious annotations")
    parser.add_argument("--fillpairs", action='store_true',
                        help="add missing annotations for mitochondrial genes")
    parser.add_argument("--standardname", action='store_true',
                        help="standardise gene names")
    parser.add_argument("--removegeneious", action='store_true',
                        help="remove any geneious editing history or other annotations")
    parser.add_argument("--dontresetsource", action = "store_true",
                        help="dont reset any source annotations to the entire sequence length")

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    # if args.filepath is not "what you want"
    #     parser.error(f"{args.filpath} is not what I want!")

    # If the arguments are all OK, output them
    return args

# Main


if __name__ == "__main__":

    # Get options

    args = getcliargs()
    #args = getcliargs("--removeduplicates".split(' '))
    #args = getcliargs("--fixhead -d INV --addanticodons --removeduplicates --fillpairs --removegeneious".split(' '))

    # Do text fixes if fixing headers
    if args.fixhead:
        furtherwork = args.fillpairs or args.standardname or args.addanticodons or \
                      args.fixlengths or args.removeduplicates or args.removegeneious
        fixout = "temp.gb" if furtherwork else sys.stdout
        fixhead(sys.stdin, fixout, args)
        #fixhead("/home/thomas/work/iBioGen_postdoc/MMGdatabase/gbmaster_2023-08-01/BIOD03404.gb", fixout, args)
    else:
        fixout = None

    if args.fillpairs or args.standardname or args.addanticodons or args.removeduplicates or \
       args.removegeneious:

        splitannotations = set()
        unrecnames = set()
        couldntduplicate = {}
        noanticodons = {}

        nameconvert, annotypes, namevariants = loadnamevariants()
        annotypes.add('gene')

        gbin = SeqIO.parse(fixout if fixout else sys.stdin, "genbank")
        for seqrecord in gbin:
            #seqrecord = list(gbin)[0]
            #sys.stderr.write(f"sequence name {seqrecord.name}\n")
            #sys.stderr.flush()
            
            # Extract all features
            #[print(f) for f in seqrecord.features[0:5]]
            
            if not args.dontresetsource:
                for i, feat in enumerate(seqrecord.features):
                    #i, feat = 0, features[0]
                    if feat.type == 'source':
                        seqrecord.features[i].location = Bio.SeqFeature.SimpleLocation(0, len(seqrecord.seq))

            if args.removegeneious:
                for feat in seqrecord.features[:]:
                    gpr = ['geneious' in qi.lower() for q in feat.qualifiers.values() for qi in q]
                    if feat.type == 'misc_feature' and any(gpr):
                        seqrecord.features.remove(feat)
                    else:
                        spr = [re.match("^[ATCGU ]+$", qi.upper()) for q in feat.qualifiers.values() for qi in q]
                        for k, gp, sp in zip(list(feat.qualifiers.keys()), gpr, spr):
                            if gp or sp or k == 'modified_by':
                                del feat.qualifiers[k]
            
            if args.removeduplicates:
                seqrecord.features = remove_duplicates(seqrecord.features)
            
            # Check for split annotation
            if args.addanticodons or args.fillpairs or args.standardname:
                locclassname = [f.location.__class__.__name__ for f in seqrecord.features]
                if any([n != "SimpleLocation" and n != "FeatureLocation" for n in locclassname]):
                    splitannotations.add(seqrecord.name)
                    sys.stdout.write(seqrecord.format('genbank'))
                    continue

            # Rename tRNAs if needed
            if args.addanticodons:
                trnas, otherfeats = sortextractfeats(seqrecord, ['tRNA', 'gene'])
                correctedtrnas, nac = add_anticodon(trnas, nameconvert)
                seqrecord.features = correctedtrnas + otherfeats
                if nac > 0:
                    noanticodons[seqrecord.name] = nac

            if args.fillpairs or args.standardname:
                # Sort features for renaming and/or duplication
                focalfeats, donefeats = sortextractfeats(seqrecord, annotypes)

                # Work through putative pairs of focal feats
                for featlis in focalfeats.values():
                    #featlis = list(focalfeats.values())[0]

                    names = [get_feat_name(f) for f in featlis]
                    types = sorted(list({f.type for f in featlis}))

                    # Iterate through the feats, making duplicates where needed
                    for i, (name, feat) in enumerate(zip(names, featlis)):
                        # i, (name, feat) = list(enumerate(zip(names, featlis)))[0]
                        # Get the standard name and type
                        stdname, ftype = [None, None]
                        if name in nameconvert:
                            stdname = nameconvert[name]
                            ftype = namevariants[stdname]['type']
                        elif name != 'unknown':
                            unrecnames.add(name)

                        # Rename if renaming and we can retrieve a standard name
                        if args.standardname and stdname:
                            feat = set_feat_name(feat, stdname)

                        # Output the renamed ones if only renaming
                        if not args.fillpairs:
                            donefeats.append(feat)
                            continue

                        # Check to see if a an appropriate duplicate is already present
                            # Are there others with this name?
                        othersamename = [j for j, n in enumerate(names) if name == n]
                        if len(othersamename) > 1:
                            # What types are they?
                            othertypes = [featlis[j].type for j in othersamename if j != i]
                            if (ftype and feat.type == 'gene' and ftype in othertypes) or (
                                    feat.type != 'gene' and 'gene' in othertypes):
                                # If the current feat is gene, with a recognised type that is
                                # present in the other features, or the feature is not a gene and
                                # gene is in the other features, sorted - just add this feature to
                                # the output
                                donefeats.append(feat)
                                continue

                        # If it's not a gene or it is a gene and we know the type to make
                        if feat.type != 'gene' or ftype:
                            newfeat = copy(feat)
                            newfeat.type = 'gene' if feat.type != 'gene' else ftype
                            donefeats.extend([feat, newfeat])

                        # Otherwise, nothing we can do because we don't know the type to make
                        else:
                            if seqrecord.name in couldntduplicate:
                                couldntduplicate[seqrecord.name] += 1
                            else:
                                couldntduplicate[seqrecord.name] = 1
                            donefeats.append(feat)
                
                seqrecord.features = donefeats

            # Output
            sys.stdout.write(seqrecord.format('genbank'))
        
        # Report any issues
        tw = _textwrap.TextWrapper(width=80, initial_indent='\t', subsequent_indent='\t')

        if len(splitannotations) > 0:
            sys.stderr.write(f"Warning: the following genbank entries had one or more annotations "
                             f"that are split over multiple locations. These should be rectified. "
                             f"Output sequences for these entries will not have anticodons added, "
                             f"pairs filled and/or names standardised.\n"
                             f"{tw.fill(', '.join(splitannotations))}\n")
        if len(noanticodons) > 0:
            nal = [f"{k}\u00A0({n})" for k, n in noanticodons.items()]
            sys.stderr.write(f"\nWarning: the following {len(noanticodons)} genbank entries "
                             f"had one or more tRNA annotations that were missing anticodon tags "
                             f"and so no anticodon could be added.\n{tw.fill(', '.join(nal))}\n")
        if len(unrecnames) > 0:
            sys.stderr.write(f"\nWarning: the following {len(unrecnames)} annotation names could "
                             f"not be recognised.\n{tw.fill(', '.join(unrecnames))}\n")
        if len(couldntduplicate) > 0:
            cdl = [f"{k}\u00A0({n})" for k, n in couldntduplicate.items()]
            sys.stderr.write(f"\nWarning: the following {len(couldntduplicate)} genbank entries "
                             f"had one or more annotations that could not be duplicated because "
                             f"the name and/or type could not be recognised.\n"
                             f"{tw.fill(', '.join(cdl))}\nNote that this may not be an issue, some "
                             f"GenBank sequences have both annotations but don't name them the "
                             f"same\n\n")
    elif args.fixlengths:
        for seqrecord in SeqIO.parse(fixout if fixout else sys.stdin, "genbank"):
            sys.stdout.write(seqrecord.format('genbank'))

    if fixout == 'temp.gb':
        os.remove('temp.gb')

exit()
