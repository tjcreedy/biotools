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


def getcliargs(arglist=None):

    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
    Standalone tool for fixing various issues in genbank files, supplied on STDIN
    |n
    To fix issues so the file can be opened with BioPython, run with --fixtext. The header will be
    checked for correct spacing, and missing division codes will be filled in with UNK by default, 
    unless a three-letter code is supplied to -d/--defaultdivision. Any incomplete annotation tags
    will be fixed as far as possible.
    |n
    On mitogenomes, fill missing annotations with --fillpairs. This will check through 
    the annotations to make sure that all 'CDS', 'tRNA' or 'rRNA' annotations have an identical 
    'gene' annotation, and vice-versa. All other types of annotation will be ignored. 
    |n 
    Standardise gene names with --standardname. This will check each gene name and replace it with
    the standard name used in https://github.com/tjcreedy/constants.   
    """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("--fixhead", action='store_true',
                        help="fix issues with headers so they work in BioPython")
    parser.add_argument("-d", "--defaultdivision", default="UNK", metavar="DIV",
                        help="a three-letter code to be used as the default devision where the "
                             "division code is missing")
    parser.add_argument("-g", "--deletegeneious", action='store_true',
                        help="try to remove automatic Geneious annotations")
    parser.add_argument("--fillpairs", action='store_true',
                        help="add missing annotations for mitochondrial genes")
    parser.add_argument("--standardname", action='store_true',
                        help="standardise gene names")

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

    # Do text fixes if fixing headers
    if args.fixhead:
        fixout = "temp.gb" if args.fillpairs or args.standardname else sys.stdout
        fixhead(sys.stdin, fixout, args)
    else:
        fixout = None

    if args.fillpairs or args.standardname:

        unrecnames = set()
        couldntduplicate = {}

        nameconvert, annotypes, namevariants = loadnamevariants()
        annotypes.add('gene')

        gbin = SeqIO.parse(fixout if fixout else sys.stdin, "genbank")

        for seqrecord in gbin:
            #seqrecord = next(gbin)

            # Find and sort features
            focalfeats, donefeats = sortextractfeats(seqrecord, annotypes)

            # Work through putative pairs of focal feats
            for featlis in focalfeats.values():
                #featlis = list(focalfeats.values())[0]

                names = [get_feat_name(f) for f in featlis]
                types = sorted(list({f.type for f in featlis}))

                # Iterate through the feats, making duplicates where needed
                for i, (name, feat) in enumerate(zip(names, featlis)):
                    #i, (name, feat) = list(enumerate(zip(names, featlis)))[0]
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
                            # If the current feat is gene, with a recognised type that is present
                            # in the other features, or the feature is not a gene and gene is in
                            # the other features, sorted - just add this feature to the output
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

            # Replace the features with the new list and output
            seqrecord.features = donefeats
            sys.stdout.write(seqrecord.format('genbank'))

        # Report any issues
        tw = _textwrap.TextWrapper(width=80, initial_indent='\t', subsequent_indent='\t')
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

        os.remove('temp.gb')

exit()
