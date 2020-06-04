#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Translate nucleotide sequences in a fasta to amino acids"""

# Imports
from Bio import Seq, SeqIO, AlignIO
import argparse
import sys
import re

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for translating nucleotide sequences in a multifasta, supplied on STDIN. All sequences are translated in the forward direction, using the same reading frame and translation table. Results are written to STDOUT")

parser.add_argument("table", help = "translation table number, required", choices = range(1,33), type = int, metavar = "TABLE")
parser.add_argument("-r","--reading_frame",help = "reading frame, default frame 1", type = int, choices = [1,2,3], default = 1)

# Function definitions

def padseq(seq, aligned = False):
    remainder = len(seq) % 3
    add = ''
    if remainder != 0:
        if aligned:
            add = '-'
        else:
            add = 'A' if seq[-remainder:] in ['T', 'TA'] else 'N'
    return(seq + add*(3-remainder))

# Main

if __name__ == "__main__":
    
    # Get options
    #args = parser.parse_args(['5'])
    args = parser.parse_args()
    
    # Read nucleotides
    # source = "/home/thomas/Documents/NHM_postdoc/MMGdatabase/0_NT_profiles_APVrev/COX1.fa"
    source = sys.stdin
    
    try: 
        nuc_records = AlignIO.read(source, "fasta")
        aligned = True
    except ValueError:
        nuc_records = SeqIO.parse(source, "fasta")
        aligned = False
    
    # Translate nucleotides
    aa_records = list()
    partialfail = dict()
    for nuc_rec in nuc_records:
        #nuc_rec = list(nuc_records)[0]
        aa_rec = nuc_rec
        # Set the correct nucleotide sequence to translate
        new_seq = nuc_rec.seq[(args.reading_frame-1):]
        # Check to see if the nucleotide sequence is formed of complete codons,
        # padding if not
        new_seq = padseq(new_seq, aligned)
        # Check for incomplete codons 
        if aligned:
            # Split into codons
            codons = [str(new_seq[i:i+3]) for i in range(0, len(new_seq), 3)]
            # Find any partial stops
            partstops = [i for i, c in enumerate(codons) if c in ['T--', 'TA-']]
            # Ensure the last partial stop is at the end of the sequence, then
            # correct it
            partstop = None
            if(len(partstops) > 0 
               and all(c == '---' for c in codons[partstops[-1] + 1:])):
                partstop = partstops[-1]
                codons[partstop] = padseq(re.sub('-', '', codons[partstop]))
            # Check for any other partials
            partial = [i for i, c in enumerate(codons) 
                       if (re.match('^-{0,2}[A-Z]{1,2}-{0,2}$', c) 
                           and i != partstop)]
            if len(partial) > 0:
                partpos = ', '.join([str(p * 3) for p in partial])
                partialfail[nuc_rec.name] = partpos
                continue
            new_seq = Seq.Seq(''.join(codons), alphabet = new_seq.alphabet)
        # Translate
        aa_rec.seq = new_seq.translate(table = args.table, gap='-')
        aa_records.append(aa_rec)
    
    if len(partialfail) > 0:
        sys.stderr.write("Warning: the following sequences had partial "
                         "internal codons and could not be translated:\n")
        for name, text in partialfail.items():
            sys.stderr.write(f"\t{name}, positions {text}\n")
    
    # Write amino acids
    SeqIO.write(aa_records, sys.stdout, "fasta")
