#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import sys
import argparse
import textwrap as _textwrap
import urllib

from collections import defaultdict

import Bio
from Bio import SeqIO

# Class definitions

# This reclasses the argparse.HelpFormatter object to have newlines in the help text for paragrahps
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

# Global variables


# Monkey Patches


Bio.SeqFeature.FeatureLocation.__hash__ = lambda self : hash(f"{self.start.__str__()}"
                                                             f"{self.end.__str__()}"
                                                             f"{self.strand}")



# Function definitions


def loadnamevariants(report=False):
    if report:
        print("The following are the gene name variants currently known.\n")
    output = {}
    fullparse = {}
    alltypes = set()
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])
        fullvariants = []
        for v in variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' '+annotype.upper()]:
                    var = v+s
                    fullvariants.append(var)
                    output[v+s] = name
        alltypes.add(annotype)
        fullparse[name] = {'type': annotype, 'variants': fullvariants}
        if report:
            fullvariants = [v.replace(' ', '\u00A0') if len(v) < 12 else v for v in fullvariants]
            fullvariants = _textwrap.fill(', '.join(fullvariants), width=80,
                                          initial_indent='\t', subsequent_indent='\t')
            print(f"Standard name = {name}, type = {annotype}, full name = {fullname}:\n"
                  f"{fullvariants}")
    if not report:
        return output, alltypes, fullparse


def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Swaps the LSU and SSU rRNA annotation locations for all entries in genbank file passed on 
        STDIN. By default, if a swap can't be made, the unchanged file will be output. 
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-s", "--strict", action='store_true',
                        help="if a swap cannot be made, do not output anything")

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # If the arguments are all OK, output them
    return args


# Start the actual script
if __name__ == "__main__":

    # Get the gene name variants
    nameconvert, annotypes, namevariants = loadnamevariants()

    # Get the arguments
    args = getcliargs() # Try to read from the command line, won't work interactively!

    # Work through records
    for seqrecord in SeqIO.parse(sys.stdin, "genbank"):
        # seqrecord = next(SeqIO.parse(source, "genbank"))

        # Sort features by location
        featdict = defaultdict(list)

        for feat in seqrecord.features:
            # feat = seqrecord.features[0]
            featdict[feat.location].append(feat)

        # Sort features into rRNAs and others
        rRNAs, outfeats = {}, []
        for loc, featlis in featdict.items():
            if any(feat.type == 'rRNA' for feat in featlis):
                rRNAs[loc] = featlis
            else:
                outfeats.extend(featlis)

        # Check the number of rRNA locations
        errmsg = None
        if len(rRNAs) == 0:
            errmsg = f"{seqrecord.name} has no rRNA locations"
        elif len(rRNAs) == 1:
            errmsg = f"{seqrecord.name} only has one rRNA annotation location"
        elif len(rRNAs) > 2:
            errmsg = f"{seqrecord.name} has more than two rRNA annotation locations"

        if errmsg:
            if args.strict:
                sys.exit(f"\nError: {errmsg}\n")
            else:
                sys.stderr.write(f"Warning: {errmsg}, no changes made\n")
                for loc, featlis in rRNAs.items():
                    outfeats.extend(featlis)
                rRNAs = {}

        # Do the rRNA swap if we still have two rRNA annotation locations
        if len(rRNAs) == 2:
            for newloc, feats in zip(reversed(rRNAs.keys()), rRNAs.values()):
                for feat in feats:
                    feat.location = newloc
                    outfeats.append(feat)

        # Output the sequence record
        seqrecord.features = outfeats
        sys.stdout.write(seqrecord.format('genbank'))

exit()
