#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports

import sys
import os
import re
import argparse
import itertools

import textwrap as _textwrap

from Bio import Align, AlignIO

# Class definitions

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width,
                                                 initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Function definitions

def informative_set(col):
    return({c.name for c in col if str(c.seq) != '-'})

def preserved_columns(palign, bycodon):
    #palign = preservealign
    ncol = palign.get_alignment_length()
    bounds = []
    for start, step in zip([0, ncol-1], [1, -1]):
        #start, step = list(zip([0, ncol-1], [1, -1]))[1]
        j = start
        while set(list(palign[:,j])) == {'-'}:
            j += step
        bounds.append(j)
    if bycodon:
        bounds = [bounds[0] - bounds[0] % 3, bounds[1] + 2 - bounds[1] % 3]
    return(bounds)

def drop_empty_rows(align):
    keepalign = []
    for a in align:
        if set(list(str(align[0].seq))) != {'-'}:
            keepalign.append(a)
    return(Align.MultipleSeqAlignment(keepalign))

def within_any(query, subjects):
    #query, subjects = infset, [r[2] for r in results]
    # s = subjects[0]
    return(any(len(query.intersection(s)) == len(query) for s in subjects))

def process_alignment(name, path, preserve, bycodon = False):
    # bycodon = args.retaincodons
    align = AlignIO.read(path, 'fasta')
    align = drop_empty_rows(align)
    preskeep = []
    preservealign = []
    if len(preserve) > 0:
        preservealign = [a for a in align if a.name in preserve]
    if len(preservealign) > 0:
        preservealign = Align.MultipleSeqAlignment(preservealign)
        preserved = preserved_columns(preservealign, bycodon)
        preskeep = range(preserved[0], preserved[1] + 1)
        checkalign = [a for a in align if a.name not in preserve]
        checkalign = Align.MultipleSeqAlignment(checkalign)
        processorder = itertools.chain(preskeep,
                                       range(0, preserved[0]),
                                       range(preserved[1] + 1,
                                             align.get_alignment_length()))
    else:
        checkalign = align
        processorder = range(align.get_alignment_length())
    del align
    processorder = list(processorder)
    results = []
    allrows = set()
    step = 3 if bycodon else 1
    for oi in range(0, len(processorder), step):
        #oi = list(range(0, len(processorder), step))[0]
        i = processorder[oi]
        infset = {s for j in range(i, i + step) 
                    for s in informative_set(checkalign[:,j:j+1])}
        if len(results) == 0:
            results.append([name, i, infset])
            allrows = infset
            if len(allrows) == len(checkalign):
                break
            continue
#        elif within_any(infset, [r[2] for r in results]):
#            continue
        elif len(allrows.union(infset)) == len(allrows):
            continue
        else:
            results.append([name, i, infset])
            allrows = allrows.union(infset)
            if len(allrows) == len(checkalign):
                break
    #[[r[0], r[1], len(r[2])] for r in results]
    results = result_filter(results, len(checkalign))
    return(results, preskeep)

def result_filter(results, maxn):
    #maxn = len(checkalign)
    results.sort(key = lambda x: len(x[2]))
    i = 0
    while i < len(results) and len(results) > 0:
        # i = 0
        if within_any(results[i][2], [r[2] for r in results[i+1:]]):
            del results[i]
        else:
            i += 1
    outresults = [results[-1]]
    i = len(results)-2
    while len(set(n for r in outresults for n in r[2])) < maxn:
        outresults.append(results[i])
        i -= 1
    return(results)

def parse_preserve(path):
    with open(path, 'r') as fh:
        out = [l.strip() for l in fh.readlines()]
    return(out)


def getcliargs(arglist = None):
    
    parser = argparse.ArgumentParser(description="""
        description:
        |n
        This script removes columns of the alignment(s) passed to 
        -a/--alignment, retaining the minimal set of columns such that all 
        sequences share information (i.e. a non gap, non ? character) with at
        least one other sequence. If multiple alignments are passed, only those
        alignments with columns required to meet this criterion are output.
        |n
        Optionally, the informative regions of some sequences can be protected
        from removal by supplying the path of a text file to -p/--preserve. The
        text file should contain the names of the sequences to preserve, one
        per line, exactly matching sequences present in at least one alignment.
        |n
        By default, column removal is arbitrary with respect to codon position.
        To retain codon position in the output alignments, i.e. retaining or 
        removing only sets of columns where each set corresponds to a complete
        codon, use the -r/--retaincodons argument. Note that no assessment of 
        reading frame is performed, codons are presumed to start at the first 
        position of each alignment and be preserved throughout with no gaps.
        |n
        Use -o/--output to specify the output. If a single alignment is
        supplied to -a/--alignment, -o/-output should be the file path to write
        the reduced alignment. If multiple alignments are supplied to 
        -a/--alignment, -o/--output should be the path to a directory into 
        which the resulting reduced alignment(s) should be written, which will
        be created if it doesn't already exist.
        """,formatter_class=MultilineFormatter)
    
    parser._optionals.title = "arguments"
    
    parser.add_argument('-a', '--alignment',
                        help = 'path to an alignment in fasta format',
                        metavar = 'path', required = True, nargs = '+',
                        type = str)
    parser.add_argument('-p', '--preserve',
                        help = 'path to a text file of sequences for which all'
                               'columns should be preserved',
                        metavar = 'path', required = False, type = str)
    parser.add_argument('-r', '--retaincodons',
                        help = 'retain codon positioning in outputs',
                        action = 'store_true')
    parser.add_argument('-o', '--output',
                        help = 'path to output file or directory',
                        required = True, metavar = 'path', type = str)
    
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    
    # Checking
    #parser.error
    
    sys.stderr.flush()
    return(args)

# Main

def main():
    
    args = getcliargs('-a genes/ND4L.fa genes/Mito_COX1_Aligned.fasta -p preserve.txt -r -o testout'.split(' '))
    
    # Create output directory if needed
    if len(args.alignment) > 1:
        os.makedirs(os.path.join(args.output), exist_ok = True)
    
    # Parse preserve list
    preserve = []
    if args.preserve:
        preserve = parse_preserve(args.preserve)
        sys.stdout.write(f"Read {len(preserve)} terminals to preserve\n")
    
    # Read and process fastas
    prescols = dict()
    results = []
    paths = dict()
    for path in args.alignment:
        #path = args.alignment[1]
        name = os.path.basename(path)
        paths[name] = path
        sys.stdout.write(f"Reducing {name}...")
        res, pc = process_alignment(name, path, preserve, args.retaincodons)
        sys.stdout.write(f"reduced to {len(res)} "
                         f"{'codon ' if args.retaincodons else ''}columns\n")
        results.extend(res)
        prescols[name] = pc
    
    # Reduce sets
    if len(args.alignment) > 1:
        sys.stdout.write(f"{len(args.alignment)} separate "
                         f"alignments produced a total of {len(results)}"
                         f"{'codon ' if args.retaincodons else ''}columns,"
                         f"now reducing these results further...\n")
        maxn = len(set(n for r in results for n in r[2]))
        results = result_filter(results, maxn)
        sys.stdout.write(f"Reduced to {len(results)} final "
                         f"{'codon ' if args.retaincodons else ''}columns "
                         f"across all alignments\n")
    
    # Write out preserved columns
    outputs = set(r[0] for r in results)
    for name in outputs:
        #name = outputs[0]
        if len(args.alignment) > 1:
            outpath = os.path.join(args.output, name)
        else:
            outpath = args.output
        pcols = set([r[1] for r in results if r[0] == name] 
                    + list(prescols[name]))
        presprin = f", including {len(prescols[name])} preserved columns,"
        sys.stdout.write(f"Writing {len(pcols)} "
                         f"{'codon ' if args.retaincodons else ''}columns"
                         f"{presprin if args.preserve else ''} ")
        step = 1
        if args.retaincodons:
            step = 3
            sys.stdout.write(f"for a total of {len(pcols) * 3} bases ")
        sys.stdout.write(f"from {name} to {outpath}\n")
        
        pcols = sorted(list(pcols))
        align = AlignIO.read(paths[name], 'fasta')
        alignkeep = align[:, pcols[0]:pcols[0]+step]
        for i in pcols[1:]:
            alignkeep += align[:,i:i+step]
        AlignIO.write(alignkeep, outpath, 'fasta')
    


if __name__ == "__main__":
    main()
    exit()

