#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Finds the frame of sequences in a multifasta"""

# Imports
from io import StringIO
from Bio import SeqIO, BiopythonWarning
import argparse
import textwrap as _textwrap
import sys
import warnings
import re
import subprocess
from copy import copy
from shutil import which


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
    chargen = (c for c in str(Seq.seq))
    gaps = 0
    while next(chargen) == '-':
        gaps += 1
    return gaps


def as_fasta(Seq):
    return f">{Seq.id}\n{str(Seq.seq)}\n"


def translate_frame(Seq, table, frame):
    # Check input
    if frame not in [1, 2, 3]:
        sys.exit("Error: frame is not 1, 2 or 3")
    # Create object
    AASeq = copy(Seq)
    # Do translation
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        Seq = Seq[(frame - 1):]
        AASeq.seq = Seq.translate(table=table).seq
    return AASeq


def quick_alignment_score(Seq):
    seqstr = str(Seq.seq)
    m = re.search(r'^-*(.*?)-*$', seqstr)
    seqstr = m.group(1)
    return seqstr.count('-') / len(seqstr.replace('-', ''))


def do_amino_acid_alignments(Seq, referencepath, frames=(1, 2, 3)):
    alignedseqs = {}
    for frame in frames:
        aaseq = translate_frame(Seq, table=5, frame=frame)
        with open("temp.fasta", 'w') as tempfasta:
            tempfasta.write(as_fasta(aaseq))
        p = subprocess.run(['mafft', '--quiet', '--6merpair', '--anysymbol',
                            '--addfragments', 'temp.fasta', referencepath],
                           capture_output=True, text=True)
        *_, aln = SeqIO.parse(StringIO(p.stdout), "fasta")
        alignedseqs[frame] = aln
    return alignedseqs


def alignment_positions(alngen, ref=0):
    """Find the position of aligned sequences in a generator relative to the sequence at the index
    supplied. Sequences prior to the supplied index are ignored"""
    # alngen, ref = SeqIO.parse(args.alignment, "fasta"), 0

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
    refstart = count_start_gaps(refseq)

    # Assess positions of subsequent sequences
    positions = dict()
    for seq in alngen:
        # seq = next(alngen)
        # Get starting position of the sequence
        seqstart = count_start_gaps(seq)
        # If the sequence starts earlier than the reference, find the number of bases earlier
        if seqstart < refstart:
            positions[seq.id] = seqstart - refstart
        # If the sequence starts later than the reference, find the number of gaps in the reference
        # up to and including this position, and subtract from the starting position
        else:
            positions[seq.id] = seqstart - str(refseq.seq[:seqstart]).count('-')

    return positions


def stopcount(seqrecord, table, frame=(1, 2, 3)):
    # Check input types
    run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
    # Run counting
    counts = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        for i in run_frame:
            seq = seqrecord.seq[(i - 1):]
            counts.append(seq.translate(table=table).count("*"))
    # Return string or list depending on length
    if (len(counts) > 1):
        return counts
    else:
        return [counts[0]]


def getcliargs(arglist=None):
    parser = argparse.ArgumentParser(description="""
        Standalone tool for finding the reading frame of each sequence in a multifasta, supplied on 
        STDIN. The sequences should be partial or complete protein coding genes. Results are 
        written to STDOUT: by default, this is the same multifasta with the frame appended to the 
        sequence headers. Alternatively, specifying -r/--report will instead write to STDOUT a tab-
        separated tabular report where the first column is headers and the subsequent column(s) the 
        found frame(s), depending on the methods selected. Reading frame can be found using three 
        methods.
        |n
        Alignment-based assessment uses an alignment of the target sequences against a reference 
        sequence or sequences to find the frame. To use this, first perform an alignment of the 
        sequences in the multifasta against one or more reference sequences that start in frame. 
        Then supply the path to this alignment to -a/--alignment. Frame will be assessed relative 
        to the start of the first sequence in the alignment, so this should be a reference 
        sequence; any subsequent sequences with headers not in the input multifasta will be assumed 
        to also be reference sequences and will be ignored.
        |n
        Stop-based assessment will be run if the -s/--stop argument is invoked. This method 
        translates the sequences read on STDIN (according to the translation table supplied to 
        -t/--table), in all three possible frames. The frame with the fewest internal stop codons 
        is selected. This method is slower and may be less accurate. 
        |n
        Reference-based assessment also translates the sequences in all three possible frames (also 
        requiring -t/--table), but then aligns each of the three translations to an amino acid 
        reference alignment or sequence supplied to -r/--reference. The frame producing the best 
        alignment (measured by the number of internal gaps created in the translation) will be 
        selected. This method is the slowest but probably the most accurate. To improve speed, this 
        method can be applied to only those input sequences that have more than one frame producing 
        the fewest stop codons, i.e. by filtering using stop-based assessment first. This is done 
        by invoking -u/--speedup and -s/--stop alongside -r/--reference. Note in this case it would 
        be somewhat redundant to give -s/--stop before -r/--reference if outputting a fasta (see 
        below).
        |n
        If arguments for more than one method (i.e. more than one of -a/--alignment, -s/--stop 
        and/or -r/--reference) are supplied, multiple methods will be used. In the case of 
        disagreements, a warning will be printed to STDERR and the precedence of the frame selected 
        will be determined by the order of the arguments supplied. If more than one method is 
        selected and -p/--report is specified, the output table will report both frames, in the 
        order of the arguments supplied; no warnings will be printed if disagreements occur.
        """, formatter_class=MultilineFormatter)

    parser.add_argument("-p", "--report", action='store_true',
                        help="write a tabular frame report to STDOUT rather than a labelled "
                             "multifasta")
    parser.add_argument("-a", "--alignment", type=str, metavar="PATH", default=argparse.SUPPRESS,
                        help="path to a nucleotide alignment comprising references and input "
                             "sequences if alignment-based assessment is required")
    parser.add_argument("-s", "--stop", action='store_true', default=argparse.SUPPRESS,
                        help="do stop-based assessment")
    parser.add_argument("-r", "--reference", type=str, metavar="PATH", default=argparse.SUPPRESS,
                        help="path to an alignment of or single amino acid reference sequence(s) "
                             "if reference-based assessment is required")
    parser.add_argument("-t", "--table", type=int, metavar="N",
                        help="translation table number if translation-based assessment and/or "
                             "reference-based assessment is required",
                        choices=range(1, 33))
    parser.add_argument("-u", "--speedup", action='store_true',
                        help="speed up reference-based assessment by filtering by stop-based "
                             "assessment first")

    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    order = [a for a in list(vars(args).keys()) if a not in ["report", 'table', 'speedup']]

    if not ('alignment' in order or 'stop' in order or 'reference' in order):
        parser.error("at least one of -a/--alignment, -s/--stop or -r/--reference is "
                     "required")
    if ('stop' in order or 'reference' in order) and not args.table:
        parser.error("-t/--table is required if using -s/--stop or -r/--reference")
    if args.speedup:
        if 'reference' not in order:
            parser.error("-s/--speedup is only relevant if using -r/--reference")
        if 'stop' not in order:
            parser.error("-s/--speedup requires -s/--stop")
        if order.index('stop') < order.index('reference') and not args.report:
            sys.stderr.write(f"Warning: -s/--stop before -r/--reference while using -s/--speedup "
                             f"will give precedence to stop-based assessment results over "
                             f"reference-based assessment results, defeating the purpose of using "
                             f"-s/--stop to pre-filter for -r/--reference\n")

    if 'reference' in order:
        if which('mafft') is None:
            sys.exit("Error: cannot locate mafft on the command line.")
    return args, order


if __name__ == "__main__":

    # Get options
    args, order = getcliargs()

    # Read nucleotides
    nuc_records = SeqIO.parse(sys.stdin, "fasta")

    # Get frame from alignment, if supplied
    alnresults = None
    if 'alignment' in order:
        positions = alignment_positions(SeqIO.parse(args.alignment, "fasta"))
        alnresults = {k: p % 3 + 1 for k, p in positions.items()}

    # Work through input sequences
    seqn = 0
    for seq in nuc_records:
        # seq = next(nuc_records)
        # Report
        seqn += 1
        sys.stderr.write(f"finding frame of sequence {seqn}\r")
        sys.stderr.flush()
        # Set up output
        frames = {}

        # Check that this sequence is present in alignment
        if 'alignment' in order:
            try:
                frames['alignment'] = {alnresults[seq.id]}
            except KeyError:
                sys.exit(f"Error: sequence {seq.id} not in {args.alignment}")

        searchframes = (1, 2, 3)

        # Get frame from translation, if requested
        if 'stop' in order:
            stops = stopcount(seq, args.table)
            frames['stop'] = {f for f, s in zip(searchframes, stops) if s == min(stops)}

        # Get frame from reference
        if 'reference' in order:
            if 'stop' in order and args.speedup:
                searchframes = frames['stop']
            if len(searchframes) > 1:
                framealns = do_amino_acid_alignments(seq, args.reference, searchframes)
                scores = [quick_alignment_score(a) for k, a in framealns.items()]
                frames['reference'] = {f for f, s in zip(searchframes, scores) if s == min(scores)}
            else:
                frames['reference'] = set(searchframes)

        # Collate and order results
        outframes = [sorted(list(frames[o]))[0] for o in order]

        # Output
        if args.report:
            print(f"{seq.id}\t", '\t'.join(str(f) for f in outframes))
        else:
            if len(set.intersection(*frames.values())) == 0:
                msg = ', '.join([f"{a} frame(s) = {v}" for a, v in frames.items()])
                sys.stderr.write(f"Warning: varying frame results for {seq.id} {msg}\n")
            print(f">{seq.id};frame={outframes[0]}\n{str(seq.seq)}")
