#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports

import sys
import argparse
import itertools

import textwrap as _textwrap

import regex
from Bio import Seq, SeqIO

# Global variables
code = {'A': 'A', 'C': 'C', 'T': 'T', 'G': 'G', 'Y': '[CT]', 'R': '[AG]',
        'W': '[AT]', 'S': '[GC]', 'K': '[TG]', 'M': '[CA]', 'D': '[AGT]',
        'V': '[ACG]', 'H': '[ACT]', 'B': '[CGT]', 'N': '[ATCG]'}


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

def required_multiple(multiple):
    class RequiredMultiple(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not len(values) % multiple == 0:
                msg = 'argument "{f}" requires a multiple of {multiple} values'
                msg = msg.format(f=self.dest, multiple=multiple)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredMultiple

def parse_primers(forward, reverse, ferror, rerror):
    #forward, reverse, ferror, rerror = args.forward, args.reverse, args.forwarderror, args.reverseerror
    
    # Add reverse complements
    primers = [forward, str(Seq.Seq(reverse).reverse_complement())]
    
    # Expand ambiguities
    primers = [''.join(code[l] for l in p) for p in primers]
    
    # Turn to regex
    primers = [regex.compile(p + e) for p, e in zip(primers, [ferror, rerror])]
    
    return(primers)

def illegal_letters(primer):
    illegal = [p for p in primer if p not in code.keys()]
    return(None if len(illegal) == 0 else illegal)


def getcliargs(arglist = None):
    
    parser = argparse.ArgumentParser(description="""
        description:
        |n
        Text
        |n
        Text
        """,formatter_class=MultilineFormatter)
    
    parser._optionals.title = "arguments"
    
    parser.add_argument('-f', '--forward',
                        help = 'forward primer sequence, 5\'-> 3\'',
                        type = str, required = True)
    parser.add_argument('-r', '--reverse',
                        help = 'reverse primer sequence, 3\'->5\'',
                        type = str, required = True)
    parser.add_argument('-i', '--input',
                        help = 'fasta file of templates',
                        type = str, required = True)
    parser.add_argument('-c', '--complete',
                        help = 'output file for complete amplicons',
                        type = str, required = True)
    parser.add_argument('-p', '--partial',
                        help = 'output file for partial amplicons, if desired',
                        type = str)
    parser.add_argument('-ef', '--forwarderror',
                        help = 'errors permitted in forward primer, in '
                               'fuzzy regex string',
                        type = str, default = '{e<=1}')
    parser.add_argument('-er', '--reverseerror',
                        help = 'errors permitted in reverse primer, in '
                               'fuzzy regex string',
                        type = str, default = '{e<=1}')
    parser.add_argument('-n', '--minlength',
                        help = 'minimum length for resulting amplicons',
                        type = int, default = 0)
    parser.add_argument('-x', '--maxlength',
                        help = 'maximum length for resulting amplicons',
                        type = int, default = float('Inf'))
    parser.add_argument('-m', '--maxamplicons',
                        help = 'maximum number of amplicons to output for '
                               'each input sequence. Input sequences with '
                               'more amplicons than this value will output 0',
                        type = int, default = 1)
    
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    
    # Checking
    args.forward, args.reverse = [p.upper()
                                  for p in [args.forward, args.reverse]]
    if illegal_letters(args.forward):
        parser.error(f"primer given to -f/--forward contains illegal letters "
                      "{','.join(illegal_letters(args.forward))}\n")
    if illegal_letters(args.reverse):
        parser.error(f"primer given to -r/--revers contains illegal letters "
                      "{','.join(illegal_letters(args.reverse))}\n")
    
    sys.stderr.flush()
    return(args)

# Main

def main():
    
    # Get arguments
    args = getcliargs()
    
    # Parse primers
    primers = parse_primers(args.forward, args.reverse, 
                            args.forwarderror, args.reverseerror)
    
    # Set up input fasta
    inseq =  SeqIO.parse(args.input, 'fasta')
    
    # Set up output fasta handles
    compfa = open(args.complete, 'w')
    if args.partial:
        partfa = open(args.partial, 'w')
    
    for template in inseq:
        #template = next(inseq)
        
        sys.stdout.write(f"{template.id} {len(template)}bp. ")
        complete = []
        partial = []
        
        for d, t in enumerate((template, template.reverse_complement())):
            #d,t = 0, template
            # Set up default sites
            sites = [[None], [None]]
            sitesn = [0, 0]
            # Do matching
            for i, p in enumerate(primers):
                #i, p = 0, primers[0]
                hits = p.finditer(str(template.seq))
                if hits:
                    sites[i] = [h.span()[1-i] for h in hits]
                    sitesn[i] = len(sites[i])
            
            sitesn = [1, 1]
            sys.stdout.write(f"{'FS' if d == 0 else 'RS'}:, "
                             f"{sitesn[0]} forward hits"
                             f"{sitesn[1]} reverse hits;")
            
            # Get amplicons
            for f, r in itertools.product(*sites):
                #f, r = sites[0][0], sites[1][0]
                if f or r:
                    amp = template[f:r]
                    if d == 1: amp.seq = amp.seq.reverse_complement()
                    if args.minlength <= len(amp) <= args.maxlength:
                        if f and r:
                            complete.append(amp)
                        else:
                            partial.append(amp)
            
            # Break out of loop if find amplicons
            if len(complete) + len(partial) > 0:
                break
        
        sys.stdout.write(f"{len(complete)} complete amplicons, "
                         f"{len(partial)} partial amplicons.\n")
        # Write out amplicons
        if 0 < len(complete) + len(partial) <= args.maxamplicons:
            for c in complete:
                compfa.write(c.format('fasta'))
            if args.partial:
                for p in partial:
                    partfa.write(p.format('fasta'))
        
    # Close handles
    compfa.close()
    if args.partial:
        partfa.close()

if __name__ == "__main__":
    main()
    exit()

