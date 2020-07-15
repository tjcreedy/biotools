#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports

import sys
import argparse
import random

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
        Randomly subsample a given proportion or number of sequences from a
        fasta file. Supply the fasta on standard input, the subsample will be
        written to standard output.
        |n
        By default all sequences are read into memory and a random set with
        count equal to the value given (if value >= 1) or the rounded product 
        of the value and the total number of input sequences (if 0 < value < 1)
        is output. For large numbers of sequences, this may be slow, so
        optionally a fast but less accurate mode can be used where the value
        given is a proportion. Each sequence has a probability of being 
        retained equal to the value given resulting in approximately the 
        desired proportion of sequences output.
        """,formatter_class=MultilineFormatter)
    parser.add_argument('value', 
                        help = ('the proportion of sequences to retain if '
                                '0 < value < 1 or the number of sequences to '
                                'retain if value is an integer => 1'),
                        type = float)
    parser.add_argument('--fast', '-f',
                        help = ('If a proportional value is given, run in '
                                'fast but less accurate mode'),
                        action = 'store_true', default = False)
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    if args.value <= 0:
        parser.error('value should be greater than 0\n')
    elif args.value >= 1:
        if args.value % 1 != 0:
            parser.error('value should be an integer if greater than 1')
        if args.fast:
            parser.error('-f/--fast is not compatable with a target number, '
                         'only a target proportion')
    sys.stderr.flush()
    return(args)

# Main

def main():
    args = getcliargs()
    seqrs = SeqIO.parse(sys.stdin, 'fasta')
    if args.fast:
        for seqr in seqrs:
            if(random.random() <= args.value):
                sys.stdout.write(seqr.format('fasta'))
    else:
        seqrs = list(seqrs)
        t = len(seqr)
        n = round(args.value * t) if args.value < 1 else args.value
        if n > t:
            sys.exit(f'Error, target subset size {n} is greater than the '
                     f'number of sequences {t}\n')
        for seqr in random.sample(seqrs, n):
            sys.stdout.write(seqr.format('fasta'))

if __name__ == "__main__":
    main()
    exit()

