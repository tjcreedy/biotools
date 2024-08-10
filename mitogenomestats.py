#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import argparse
import textwrap as _textwrap
import urllib
import os
import sys
from collections import defaultdict
from statistics import median

from Bio import __version__ as bioversion
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio import SeqFeature





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

stdorder = ('ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND5', 'ND4', 'ND4L', 'ND6', 
            'CYTB', 'ND1')

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


def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


def istrunc(feat):
    s = not isinstance(feat.location.start, SeqFeature.ExactPosition)
    e = not isinstance(feat.location.end, SeqFeature.ExactPosition)
    return s, e

def get_codons(feat, parent_sequence):
    seq = feat.extract(parent_sequence)
    start, cs, stop = str(seq[:3]), 0, ''
    if istrunc(feat)[0]:
        start = ''
        if 'codon_start' in feat.qualifiers:
            cs = int(feat.qualifiers['codon_start'][0])-1
        else:
            sys.stderr.write(f"Warning: CDS feature {get_feat_name(feat)} in {seqrecord.name} has a "
                              "truncated start but no codon_start - assuming 1\n")
            cs = 1
    if not istrunc(feat)[1]:
        remainder = len(seq[cs:])%3
        stop = str(seq[-remainder:] if remainder > 0 else seq[-3:])
    return start, stop


def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
    This script parses annotated mitogenomes stored in a genbank-format flat file and returns one 
    or two tables describing various features of the mitogenome, annotations and sequence. Supply 
    a path to -m/--mitostats to generate a summary table for mitogenomes, and/or a path to 
    -g/--genestats to generate a summary table for gene annotations. Currently only protein coding 
    genes are reported.
    """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-i", "--genbank", type=str, metavar='PATH', required=True, nargs='+',
                        help="the path(s) to one or more genbank files from which gene sequences "
                             "should be extracted")
    parser.add_argument("-m", "--mitostats", type=str, metavar='PATH', required=True,
                        help="the path to which mitogenome statistics csv will be written, if "
                             "desired")
    parser.add_argument("-g", "--genestats", type=str, metavar='PATH', required=True,
                        help="the path to which gene statistics csv will be written, if "
                             "desired")

    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # Do some checking of the inputs
    if not args.mitostats and not args.genestats:
        parser.error(f"supply one or both of -m/--mitostats and -g/--genestats")

    # If the arguments are all OK, output them
    return args


# Start the actual script
if __name__ == "__main__":

    # Get the gene name variants
    nameconvert, annotypes, namevariants = loadnamevariants()

    # Get the arguments
    args = getcliargs(None)
    #args = getcliargs('-i /home/thomas/work/iBioGen_postdoc/MMGdatabase/gbmaster_2024-07-24/GBDL02016.gb -g /home/thomas/scratch/genes.csv -m /home/thomas/scratch/mito.csv'.split(' '))  # Read from a string, good for testing
    sys.stderr.write(f"Biopython version {bioversion}\n")
    # Set up for unrecognised genes
    unrecgenes = defaultdict(list)

    # Set up output handles
    if args.mitostats:
        args.mitostats = open(args.mitostats, 'w')
        args.mitostats.write("mt_id,length,ngenes,topology,GCfraction,genes,genepresence,geneorder\n")
        args.mitostats.flush()
    if args.genestats:
        args.genestats = open(args.genestats, 'w')
        args.genestats.write("mt_id,gene,GCfraction,prev_genename,prev_genedist,next_genename,"
                             "next_genedist,startcodon,starttrunc,stopcodon,stoptrunc\n")
        args.genestats.flush()


    # Work through genbank files
    for gbpath in args.genbank:
        # gbpath = args.genbank[0]
        for seqrecord in SeqIO.parse(gbpath, "genbank"):
            # gb = SeqIO.parse(gbpath, "genbank")
            # seqrecord = next(gb)
            
            # Extract names
            seqname = seqrecord.name

            # Set up genes containers
            genes = {}
            
            # Work through features
            for feat in seqrecord.features:
                # feat = seqrecord.features[3]
                # Convert feat name
                name = get_feat_name(feat)
                if name in nameconvert:
                    stdname = nameconvert[name]
                else:
                    unrecgenes[name].append(seqname)
                    continue

                if name in stdorder:
                    if name in genes and genes[name].type == 'CDS':
                        continue
                    else:
                        genes[name] = feat

            # Gene order
            genes = [(n, f) for n, f in genes.items()]
            if len(genes) > 1:
                    # Sort by location on mitogenome
                genes = sorted(genes, key=lambda i: int(i[1].location.start))
                order = [i[0] for i in genes]
                    # Find direction of genes and flip if necessary
                iorder = [stdorder.index(l) for l in order]
                diffs = [y - x for x, y in zip(iorder, iorder[1:])]
                if median(diffs) < 0:
                    genes.reverse()
                    order = [i[0] for i in genes]
                    iorder = [order.index(l) for i, l in enumerate(stdorder) if l in order]
                    diffs = [y - x for x, y in zip(iorder, iorder[1:])]
                # Reset to a common starting point
                if min(diffs) < -10:
                    breaki = diffs.index(min(diffs))+1
                    genes = genes[breaki:] + genes[:breaki]

            if args.mitostats:
                args.mitostats.write(','.join([
                    seqname,
                    str(len(seqrecord.seq)),
                    str(len(genes)),
                    seqrecord.annotations['topology'],
                    str(gc_fraction(seqrecord.seq)),
                    ''.join('X' if i in order else '-' for i in stdorder),
                    ';'.join(order)
                    ]) + '\n')

            # Gene statistics
            if args.genestats:
                for i, (n, f) in enumerate(genes):
                    # i, (n, f) = list(enumerate(genes))[0]
                    seq = f.extract(seqrecord.seq)
                    gene_codons = get_codons(f, seqrecord.seq)
                    
                    prev_genename, next_genename = '', ''
                    prev_genedist, next_genedist = '', ''
                    if not istrunc(f)[0] and i > 0:
                        bi = i -1
                        prev_genename = genes[bi][0]
                        if not istrunc(genes[bi][1])[1]:
                            be = int(genes[bi][1].location.end)
                            if int(genes[bi][1].location.start) < int(f.location.start):
                                prev_genedist = int(f.location.start) - be
                            else: 
                                prev_genedist = int(f.location.start) + len(seqrecord.seq) - be
                        else:
                            prev_genedist = ''
                    if not istrunc(f)[1] and i < len(genes)-1:
                        ni = i + 1
                        next_genename = genes[ni][0]
                        if not istrunc(genes[ni][1])[0]:
                            ns = int(genes[ni][1].location.start)
                            if int(genes[ni][1].location.end) > int(f.location.end):
                                next_genedist = ns - int(f.location.end)
                            else:
                                next_genedist = ns + len(seqrecord.seq) - int(f.location.end)
                        else:
                            next_genedist = ''

                    args.genestats.write(','.join([
                        seqname, n, str(len(seq)), str(gc_fraction(seq)),
                        prev_genename, str(prev_genedist), next_genename, str(next_genedist),
                        gene_codons[0], str(istrunc(f)[0]), gene_codons[1], str(istrunc(f)[1])
                        ]) + '\n')
                        


    # Close file handles
    if args.mitostats:
        args.mitostats.close()
    if args.genestats:
        args.genestats.close()

    # Report unrecognised genes
    if len(unrecgenes) > 0:
        sys.stderr.write(f"Warning: unrecognised genes present in the following entries.\n")
        for gene, entries in unrecgenes.items():
            sys.stderr.write(f"\t{gene} - {', '.join(sorted(list(set(entries))))}\n")

exit()
