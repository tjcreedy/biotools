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
    parser.add_argument('-a', '--padto',
                        help = 'add Ns to the start or end of partial '
                               'amplicons shorter than this value to make up '
                               'to this length',
                        type = int, default = 0)
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
    #args = getcliargs('-i /home/thomas/Downloads/MIDORI/MIDORI_UNIQUE/MIDORI_UNIQUE_20180221_COI_MOTHUR.fasta -f CCNGAYATRGCNTTYCCNCG -r TANACYTCNGGRTGNCCRAARAAYCA -c MIDORI_UNIQUE_20180221_COI_418_complete.fasta -p MIDORI_UNIQUE_20180221_COI_418_partial.fasta -n 350 -x 436 -m 1 -a 418'.split(' '))
    args = getcliargs()
    
    # Parse primers
    primers = parse_primers(args.forward, args.reverse, 
                            args.forwarderror, args.reverseerror)
    
    # Set up input fasta
    inseq =  SeqIO.parse(args.input, 'fasta')
    #indict = SeqIO.to_dict(inseq)
    
    # Set up output fasta handles
    compfa = open(args.complete, 'w')
    if args.partial:
        partfa = open(args.partial, 'w')
    
    for template in inseq:
        #template = indict['AB564645.1._1._1044']
        #template = next(inseq)
        
        sys.stdout.write(f"{template.id} {len(template)}bp. ")
        complete = []
        partial = []
        
        for d, t in enumerate((template, template.reverse_complement())):
            #d,t = 0, template
            #d,t = 1, template.reverse_complement()
            
            # Generate hit lists for each primer
            hits = [p.finditer(str(t.seq)) for p in primers]
            # Extract start/stop positions for spanned regions
            sites = [[h.span()[1-i] for h in hs] for i, hs in enumerate(hits)]
            # Count hits
            sitesn = [len(s) for s in sites]
            # Replace empty lists with None
            sites = [s if len(s) > 0 else [None] for s in sites]
            
            sys.stdout.write(f"{'FS' if d == 0 else 'RS'}: "
                             f"{sitesn[0]} forward hit(s), "
                             f"{sitesn[1]} reverse hit(s), ")
            
            # Get amplicons
            for f, r in itertools.product(*sites):
                #f, r = list(itertools.product(*sites))[0]
                lens = []
                if f or r:
                    amp = template[f:r]
                    lens.append(str(len(amp)))
                    if d == 1: amp.seq = amp.seq.reverse_complement()
                    if args.minlength <= len(amp) <= args.maxlength:
                        if f and r:
                            complete.append(amp)
                        else:
                            if 0 < args.padto and len(amp) < args.padto:
                                pad = Seq.Seq('N' * (args.padto - len(amp)))
                                if (f and d == 0) or (r and d == 1):
                                    amp.seq = amp.seq + pad
                                else:
                                    amp.seq = pad + amp.seq
                            partial.append(amp)
            
            sys.stdout.write(f"{', '.join(lens)}bp amplicon(s); ")
            # Break out of loop if found any hits (even if rejected)
            if len(lens) > 0:
                break
        
        sys.stdout.write(f"written {len(complete)} complete amplicon(s), "
                         f"{len(partial)} partial amplicon(s)\n")
        
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

