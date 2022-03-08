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


# Argument parser

def getcliargs():
    parser = argparse.ArgumentParser(description="""
        description:
        |n
        Standalone tool for filtering an alignment supplied in a multifasta on
        STDIN. Currently two options: remove any rows matching a regular 
        expression supplied to --removeseqs, e.g. to remove any sequences 
        beginning 'NT_', pass --removeseqs '^NT_'. To remove any columns that 
        are only gaps, pass --removegapcols. This is run after any sequence
        removal. The output multifasta alignment is written to STDOUT""",
                                     formatter_class=MultilineFormatter)
        
    parser._optionals.title = "arguments"
    
    parser.add_argument('--removeseqs', metavar = 'RX', 
                        help = 'remove sequences matching this regex',
                        type = str)
    parser.add_argument('--removegapcols', 
                        help = 'remove columns that are only gaps',
                        action = 'store_true')
    
    return(parser.parse_args())

# Main

def main():
    args = getcliargs()
    # Read in alignment
    aln = AlignIO.read(sys.stdin, 'fasta')
    
    if args.removeseqs:
        row_keep = []
        for a in aln:
            if not re.match(args.removeseqs, a.name):
                row_keep.append(a)
        aln = Align.MultipleSeqAlignment(row_keep)
    
    
    if args.removegapcols:
        ncol = aln.get_alignment_length()
        col_keep = []
        for i in range(ncol):
            if set(list(aln[:,i])) != {'-'}:
                col_keep.append(i)
        aln_keep = aln[:,col_keep[0]:col_keep[0]+1]
        for i in col_keep[1:]:
            aln_keep += aln[:,i:i+1]
        aln = aln_keep
    
    AlignIO.write(aln, sys.stdout, 'fasta')


if __name__ == "__main__":
    main()
    exit()

