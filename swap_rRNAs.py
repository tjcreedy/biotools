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
        STDIN 
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    # parser.add_argument("-g", "--genbank", type=str, metavar='PATH', required=True, nargs='+',
    #                     help="the path(s) to one or more genbank files")

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
        if len(rRNAs) == 0:
            sys.stderr.write(f"{seqrecord.name} has no rRNA locations, no changes made\n")
        elif len(rRNAs) == 1:
            sys.stderr.write(f"{seqrecord.name} only has one rRNA annotation location, no changes "
                             f"made\n")
            outfeats.extend(list(rRNAs.values())[0])
            rRNAs = {}
        elif len(rRNAs) > 2:
            sys.exit(f"\n\nERROR: {seqrecord.name} has more than two rRNA annotation locations, "
                     f"cannot perform a swap!\n\n")

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
