#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports

import sys
import argparse

import textwrap as _textwrap

from Bio import SeqIO

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
def getcliargs(arglist = None):
    parser = argparse.ArgumentParser(description="""
        description:
        |n
        Filter the sequences in a fasta based on their exact presence in 
        a second fasta. Supply the fasta to be filtered on standard input, 
        filtered sequences will be written to standard output.
        """,formatter_class=MultilineFormatter)
    parser.add_argument('action',
                        help = ('specify whether to \'k\'eep or \'r\'emove '
                                'sequences found in filter'),
                        type = str, choices = ['k', 'r'])
    parser.add_argument('filter',
                        help = ('path to the fasta containing the filter '
                                'sequences'),
                        type = str)
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    sys.stderr.flush()
    return(args)


# Main

def main():
    args = getcliargs()
    seqfilter = {str(s.seq) for s in SeqIO.parse(args.filter, 'fasta')}
    for seqr in SeqIO.parse(sys.stdin, 'fasta'):
        seq = str(seqr.seq)
        if ((seq in seqfilter and args.action == 'k')
            or (seq not in seqfilter and args.action == 'r')):
            sys.stdout.write(seqr.format('fasta'))

if __name__ == "__main__":
    main()
    exit()

