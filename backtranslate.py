#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Back-translate alignment of amino acids in a fasta to aligned nucleotide sequences using unaligned nucleotide sequences as a guide"""

# Imports
from Bio import SeqIO
import argparse
import sys
import re

# Global variables

parser = argparse.ArgumentParser(
    description="""
    Standalone tool for converting aligned amino acid sequences and the source unaligned 
    nucleotide sequences into the equivalent aligned nucleotide sequences. All sequences are 
    back-translated in the forward direction, using the same translation table. By default, if  
    sequences in one of the files have a ";frame=N" tag in the sequence header, where N is 1, 2 or 
    3, each sequence will be backtranslated according to this frame. Otherwise, sequences will be 
    backtranslated in frame 1. To override this, specify a frame with -r/--readingframe.
    Amino acid alignment should be supplied on STDIN, unaligned nucleotides should be supplied as 
    a file path. Unaligned and aligned entries should be in exactly the same order. Results are 
    written to STDOUT
    """)

parser.add_argument("nucpath", help="path to corresponding unaligned nucleotides",
                    metavar="NUCPATH")
parser.add_argument("table", help="translation table number, required", choices=range(1, 33),
                    type=int, metavar="TABLE")
parser.add_argument("-i", "--input", help="amino acid alignment", type=str)
parser.add_argument("-r", "--readingframe", help="reading frame override", type=int,
                    choices=[1, 2, 3])
parser.add_argument("-s", "--striptags", help="remove frame tags in the output file",
                    action='store_true')


# Function definitions

def stripframe(st):
    return re.sub(r';frame=[0-9]+(;$)?', '', st)


def count_start_gaps(Seq):
    chargen = (c for c in str(Seq.seq))
    gaps = 0
    while next(chargen) == '-':
        gaps += 1
    return gaps


def count_end_gaps(Seq):
    return count_start_gaps(Seq[::-1])



# TODO: there's a module in biopython for doing this! https://biopython.org/DIST/docs/api/Bio.codonalign-module.html

# Main

if __name__ == "__main__":

    # Get options

    args = parser.parse_args()

    # Read amino acid alignment
    if (args.input):
        aaa_records = SeqIO.parse(args.input, "fasta")
    else:
        aaa_records = SeqIO.parse(sys.stdin, "fasta")
    # Read nucleotides
    nuc_records = SeqIO.parse(args.nucpath, "fasta")

    # Compare lengths of generators
    # sys.exit("Error: number of sequences in amino acid alignment (", len(aaa_records), ") does not match number of sequences in unaligned nucleotides file (", len(nuc_records), ")")

    # Find gap positions and convert nucleotide data
    n = 0
    nuc_aln = list()
    taggedrf = None
    waste = [float('Inf')] * 2
    for aaa_rec, nuc_rec in zip(aaa_records, nuc_records):
        # aaa_rec, nuc_rec = list(zip(aaa_records, nuc_records))[14]

        # sys.stderr.write("Comparing AA %s %s with NT %s %s\n" % (aaa_rec.id, aaa_rec.seq, nuc_rec.id, nuc_rec.seq))

        n += 1

        # Run checks
        if stripframe(aaa_rec.id) != stripframe(nuc_rec.id):
            sys.exit("Error: sequence identifiers do not match for sequence " + str(
                n) + " (alignment: " + aaa_rec.id + " nucleotides: " + nuc_rec.id)

        # Get the reading frame
        if args.readingframe is None:
            rf = []
            # Check tags from both input files and retrieve RF(s) if present
            for r in [aaa_rec, nuc_rec]:
                rftag = re.search("(?<=;frame=)([0-9]+)", r.id)
                if rftag:
                    rf.append(int(rftag.group(1)))
            rf = set(rf)
            if len(rf) > 1:
                sys.exit(
                    f"Error: both input fastas have frame tags, but they do not agree for "
                    f"sequence {n}: AA={aaa_rec.id}, NT={nuc_rec.id}")
            elif len(rf) == 0:
                rf = None
            else:
                rf = list(rf)[0]

            mixerrmsg = "Error: some, but not all sequences have \";frame=\" tag, either " \
                        "correct this or run with -r/--readingframe to ignore tags and run with " \
                        "a universal frame"
            if rf:
                if taggedrf is not None and not taggedrf:
                    sys.exit(mixerrmsg)
                taggedrf = True
                if rf not in (1, 2, 3):
                    sys.exit(f"Error: frame tag for sequence {nuc_rec.id} not recognised")
                if args.striptags:
                    for a in ['id', 'name', 'description']:
                        setattr(nuc_rec, a, stripframe(getattr(nuc_rec, a)))
            elif taggedrf:
                sys.exit(mixerrmsg)
            else:
                rf = 1
        else:
            rf = args.readingframe

        # Convert to 0-indexed reading frame
        rf0 = rf - 1

        # Add gaps for incomplete ends
        remainder = len(nuc_rec.seq[(rf0):]) % 3
        nuc_rec.seq = nuc_rec.seq + '-' * (3 - remainder) if (remainder != 0) else nuc_rec.seq

        # Find gap positions
        aa_gaps = [i for i, b in enumerate(aaa_rec.seq) if b == "-"]

        # Extract leading gaps
        starti = 0
        for starti in range(len(aa_gaps)):
            if not (aa_gaps[starti] == 0 or (aa_gaps[starti] - aa_gaps[starti-1] == 1)):
                break
        aa_gaps = aa_gaps[starti:]

        # Convert internal gaps to nucleotide positions, offset by reading frame
        nuc_gaps = [a * 3 + rf0 for a in aa_gaps]

        # Add leading gaps
        nuc_rec.seq = '-' * (3 * starti) + nuc_rec.seq

        # Insert gaps
        for g in nuc_gaps:
            nuc_rec.seq = nuc_rec.seq[:g] + "---" + nuc_rec.seq[g:]

        # Add trailing gaps if frame requires it and it's less than 3
        lendiff = len(str(aaa_rec.seq)) * 3 - (len(nuc_rec.seq) + rf0)
        if (lendiff != 0):
            if (lendiff < 3):
                nuc_rec.seq = nuc_rec.seq + ("-" * lendiff)
            else:
                sys.exit(
                    "Error: number of bases / amino acids do not correspond in sequence number " + str(
                        n) + " " + aaa_rec.id + " (alignment)")

        # Pad ahead of sequence to ensure in frame if we have/get reading frames other than 1
        nuc_rec.seq = "-" * (3 - rf0) + nuc_rec.seq

        # Count number of gaps at start and end
        padgaps = [f(nuc_rec) for f in [count_start_gaps, count_end_gaps]]
        for i in range(2):
            waste[i] = min([padgaps[i], waste[i]])

        nuc_aln.append(nuc_rec)

    # Remove any unecessary leading or trailing gaps, if necessary
    waste = [int(w/3)*3 for w in waste]
    if any(w > 0 for w in waste):
        for i in range(len(nuc_aln)):
            nuc_aln[i].seq = nuc_aln[i].seq[waste[0]:len(nuc_aln[i].seq)-waste[1]]

    # Write nucleotides
    SeqIO.write(nuc_aln, sys.stdout, "fasta")
