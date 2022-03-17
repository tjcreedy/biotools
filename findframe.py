#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Translate nucleotide sequences in a fasta to amino acids"""

# Imports
from Bio import Seq, SeqIO, AlignIO, BiopythonWarning
import argparse
import sys
import warnings
import re


# Class definitions

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



# Function definitions

def count_start_gaps(Seq):
    chargen = (c for c in str(Seq))
    gaps = 0
    while next(chargen) == '-':
        gaps += 1
    return gaps

def do_amino_acid_alignments(seq, referencepath):

    # Create temporary directory





def alignment_positions(alngen, ref = 0):
    """Find the position of aligned sequences in a generator relative to the sequence at the index
    supplied. Sequences prior to the supplied index are ignored"""
    #alngen, ref = SeqIO.parse(args.alignment, "fasta"), 0

    # Skip any sequences up to the reference and get the reference sequence
    i = 0
    refseq = None
    while i <= ref:
        i += 1
        try:
            refseq = next(alngen)
        except StopIteration:
            sys.exit("Error: ref is too high, not enough values in alngen")

    # Get the starting position of the reference
    refstart = count_start_gaps(refseq.seq)

    # Assess positions of subsequent sequences
    positions = dict()
    for seq in alngen:
        #seq = next(alngen)
        # Get starting position of the sequence
        seqstart = count_start_gaps(seq.seq)
        # If the sequence starts earlier than the reference, find the number of bases earlier
        if seqstart < refstart:
            positions[seq.id] = seqstart - refstart
        # If the sequence starts later than the reference, find the number of gaps in the reference
        # up to and including this position, and subtract from the starting position
        else:
            positions[seq.id] = seqstart - str(refseq.seq[:seqstart]).count('-')

    return positions


def stopcount(seqrecord, table, frame = (1 ,2, 3)):
    # Check input types
    run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
    # Run counting
    counts = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        for i in run_frame:
            seq = seqrecord.seq[(i-1):]
            counts.append(seq.translate(table = table).count("*"))
    # Return string or list depending on length
    if(len(counts) > 1):
        return counts
    else:
        return [counts[0]]

def getcliargs(arglist=None):
    parser = argparse.ArgumentParser(description="""
        Standalone tool for finding the reading frame of each sequence in a multifasta, supplied on 
        STDIN. The sequences are assumed to be partial or complete protein coding genes. Results 
        are written to STDOUT: by default, this is the same multifasta with the frame appended to 
        the sequence headers. Alternatively, specifying -r/--report will instead write to STDOUT a 
        tab-separated tabular report where the first column is headers and the subsequent column(s) 
        the found frame(s), depending on the methods selected.
        |n
        Reading frame can be found using two methods. Reference-based assessment uses an alignment 
        of the target sequences against a reference sequence or sequences to find the frame. To use 
        this, perform an alignment of the sequences in the multifasta against one or more reference 
        sequences that start in frame and supply the path to this alignment to -a/--alignment. 
        Frame will be assessed relative to the start of the first sequence in the alignment, so 
        this should be a reference sequence; any subsequent sequences with headers not in the input 
        multifasta will be assumed to also be reference sequences and will be ignored.
        |n
        Translation-based assessment translates the sequences read on STDIN according to the 
        translation table supplied to -t/--table, in all three possible frames. The frame with the 
        fewest internal stop codons is selected. This method is slower and may be less accurate. 
        |n
        If both -a/--alignment and -t/--table are supplied, both methods will be used. In the case 
        of disagreements, a warning will be printed to STDERR and the precedence of the frame 
        selected will be determined by the order of the arguments supplied. If both methods are 
        selected and -r/--report is specified, the output table will report both frames, in the 
        order of the arguments supplied; no warnings will be printed if disagreements occur.
        """, formatter_class=MultilineFormatter)

    parser.add_argument("-r", "--report", action='store_true',
                        help="write a tabular frame report to STDOUT rather than a labelled multifasta")
    parser.add_argument("-t", "--table", type=int, metavar="TABLE", default=argparse.SUPPRESS,
                        help="translation table number if translation-based assessment is required",
                        choices=range(1, 33))
    parser.add_argument("-a", "--alignment", type=str, metavar="PATH", default=argparse.SUPPRESS,
                        help="path to an alignment if reference-based assessment is required")

    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    order = [a for a in list(vars(args).keys()) if a != "report"]

    if not (args.table or args.alignment):
        parser.error("at least one of -t/--table and -a/--alignment is required")

    return args, order


if __name__ == "__main__":

    # Get options
    #args, order = getcliargs('-t 5'.split(' '))
    args, order = getcliargs()

    # Read nucleotides
    # nuc_records = SeqIO.parse("test/sequence.fasta", 'fasta')
    nuc_records = SeqIO.parse(sys.stdin, "fasta")

    # Get frame from alignment, if supplied
    refframe = None
    if 'alignment' in order:
        positions = alignment_positions(SeqIO.parse(args.alignment, "fasta"))
        refframe = {k: p%3 + 1 for k, p in positions.items()}

    # Work through input sequences
    for seq in nuc_records:
        #seq = next(nuc_records)

        # Check that this sequence is present in alignment
        if 'alignment' in order:
            if seq.id not in refframe.keys():
                sys.exit(f"Error: sequence {seq.id} not in {args.alignment}")

        # Get frame from translation, if requested
        transframe = None
        if 'table' in order:
            stops = stopcount(seq, args.table)
            transframe = min(range(len(stops)), key=stops.__getitem__) + 1

        # Collate and order results
        frames = [transframe, refframe[seq.id] if 'alignment' in order else None]
        if order[0] == "alignment":
            frames.reverse()
        frames = [f for f in frames if f]

        # Output
        if args.report:
            print(f"{seq.id}\t", '\t'.join(str(f) for f in frames))
        else:
            if len(frames) == 2 and frames[0] != frames[1]:
                sys.stderr.write(f"Warning: {seq.id} reference frame = {str(refframe[seq.id])} but "
                                 f"translation frame = {str(transframe)}\n")
            print(f">{seq.id};frame={frames[0]}\n{str(seq.seq)}")