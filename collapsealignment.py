#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Description of what the script does"""

# Imports

import argparse
import sys
from collections import defaultdict, Counter
import textwrap as _textwrap



from Bio import Seq, SeqRecord, Align, AlignIO
from Bio.Align import AlignInfo

# Class definitions

# This reclasses the argparse.HelpFormatter object to have newlines in the help text for paragrahps
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

def pairwise_dist(a, b):
    c = []
    for ai, bi in zip(str(a).upper(), str(b).upper()):
        if not(ai == "-" or bi == "-"):
            if ai == bi or ai == "X" or bi == "X":
                c.append(1)
            else:
                c.append(0)
    if len(c) == 0:
        return None
    else:
        return sum(c)/len(c)

def float_range(mini,maxi):
    """Return function handle of an argument type function for
       ArgumentParser checking a float range: mini <= arg <= maxi
         mini - minimum acceptable argument
         maxi - maximum acceptable argument"""

    # Define the function with default arguments
    def float_range_checker(arg):
        """New Type function for argparse - a float within predefined range."""

        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < mini or f > maxi:
            raise argparse.ArgumentTypeError("must be in range [" + str(mini) + " .. " + str(maxi)+"]")
        return f

    # Return function handle to checking function
    return float_range_checker


def getcliargs(arglist=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
        Collapse multiple sequences with the same name in an alignment into one sequence.  
        |n
        If sequences with the same name overlap with each other, and the overlapping region has an 
        edit distance less than -t/--threshold, the sequence most dissimilar (edit distance) from 
        the consensus of all sequences with unique names is discarded. If the edit distance of the 
        overlapping distance is greater that -t/--threshold, the majority consensus will be used. 
        """, formatter_class=MultilineFormatter)

    # Add individual argument specifications
    parser.add_argument("-i", "--input", type=str, metavar="PATH", required=True,
                        help="path to an input alignment")
    parser.add_argument("-t", "--threshold", type=float_range(0,1), metavar="[0-1]", default=1,
                        help="minimum edit distance of overlapping regions to retain sequences")
    parser.add_argument("-o", "--output", type=str, metavar="PATH", required=True,
                        help="path to write collapsed alignment")
    # Parse the arguments from the function call if specified, otherwise from the command line
    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    # If the arguments are all OK, output them
    return args


# Start the actual script
if __name__ == "__main__":
    # Get the arguments
    args = getcliargs()  # Try to read from the command line, won't work interactively!

    # Open alignment
    aln = AlignIO.read(args.input, 'fasta')
    print(f"Read alignment of {len(aln)} sequences")

    # Create hash of alignment records by name
    records = defaultdict(list)
    for a in aln:
        records[a.name].append(a)

    # Subset alignment to only records with a single sequence
    refnames = [r for r,s in records.items() if len(s) == 1]
    refaln = Align.MultipleSeqAlignment(a for a in aln if a.name in refnames)

    print(f"Alignment has {len(refaln)} sequences with unique names")

    # Subset dict to only non-reference sequences
    records = {r: s for r, s in records.items() if len(s) > 1}

    # Generate consensus of references
    consensus = Align.AlignInfo.SummaryInfo(refaln).gap_consensus()

    # Score non-reference sequences in similarity to consensus
    scores = defaultdict(list)
    for r, s in records.items():
        print(f"Found {len(s)} sequences with name {r}")
        # r,s = list(records.items())[0]
        for seq in s:
            # seq = s[0]
            scores[r].append(pairwise_dist(consensus, seq.seq))

    # Find overlapping sequences and remove the worst ones
    for r in records.keys():
        # r = list(records.keys())[0]
        s = records[r]
        # Find the overlaps - sequences that have a pd must overlap, but want to retain
        # both sequences if the overlap is very similar
        overlaps = []
        for i in range(len(s)):
            for j in range(i+1, len(s)):
                # i,j = 0,1
                pd = pairwise_dist(s[i].seq, s[j].seq)
                if pd is not None and pd < args.threshold:
                    overlaps.append((i, j))
        print(f"{r} sequences have {len(overlaps)} overlapping regions below threshold")
        # Iteratively reject the sequence with the worst similarity to the consensus
        # until no overlaps remain
        reject = []
        while len(overlaps) > 0:
            queries = list({o for ol in overlaps for o in ol})
            qscore = [scores[r][q] for q in queries]
            rj = [q for q,qs in zip(queries, qscore) if qs == min(qscore)][0]
            reject.append(rj)
            overlaps = [o for o in overlaps if rj not in o]
        print(f"Discarding {len(reject)} sequences from {r}")
        # Drop the rejects
        s = [seq for i,seq in enumerate(s) if i not in reject]
        records[r] = s

    # Collapse the sequences
    collapsedseq = []
    for r, s in records.items():
        # r,s = list(records.items())[0]
        print(f"Collapsing {len(s)} sequences from {r}")
        raln = Align.MultipleSeqAlignment(s)
        consensus = ""
        for n in range(raln.get_alignment_length()):
            bases = Counter()
            for rec in raln:
                b = rec[n]
                if b != '-':
                    bases[b] += 1
            if len(bases) == 0:
                consensus += "-"
            elif len(bases) == 1:
                consensus += list(bases.keys())[0]
            else:
                counts = list(bases.values())
                mfreq = [b for b,c in bases.items() if c == max(counts)]
                if len(mfreq) == 1:
                    consensus += mfreq[0]
                else:
                    consensus = "N"
        collapsedseq.append(SeqRecord.SeqRecord(Seq.Seq(consensus), id=r, description=''))

    for c in collapsedseq:
        refaln.append(c)
    print(f"Finished: writing alignment with {len(refaln)} sequences")
    AlignIO.write(refaln, args.output, 'fasta')



