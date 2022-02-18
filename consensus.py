#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate consensus sequence from alignment"""

# Imports
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
import argparse
import sys

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

# Global variables

dna_to_ambig = {'A' : 'A',
         'C' : 'C',
         'G' : 'G',
         'T' : 'T',
         'AC' : 'M',
         'AG' : 'R',
         'AT' : 'W',
         'CG' : 'S',
         'CT' : 'Y',
         'GT' : 'K',
         'ACG' : 'V',
         'ACT' : 'H',
         'AGT' : 'D',
         'CGT' : 'B',
         'ACGT' : 'N'
        }
ambig_to_dna = dict()
for dna in dna_to_ambig:
    if dna in list('ATCG'):
        continue
    else:
        ambig_to_dna[dna_to_ambig[dna]] = dna


parser = argparse.ArgumentParser(description = "Standalone tool for generating a consensus sequence from an alignment in a multifasta, supplied on STDIN. Results are written to STDOUT")

parser.add_argument("-t","--threshold", help = "threshold for selecting residues", type = float_range(0, 1), default = 0.5)
parser.add_argument("-n","--name", help = "name of the output sequence", type = str, default = "consensus")

# Function definitions

threshold=.3
def mr_consensus(aln, threshold=0.7):
        """Output a majority rule consensus sequence of the alignment, allowing
        gaps. For each position, if any nucleotide frequencies exceed 
        threshold, the nucleotide, gap, or appropriate ambiguity (threshold 
        <.5) is returned; otherwise, the appropriate ambiguity code is returned
        based on all present nucleotides. If gaps are present, they will be 
        output if more frequent than all nucleotides
        """
        
        alinf = AlignInfo.SummaryInfo(aln)
        consensus = ""

        # find the length of the consensus we are creating
        con_len = alinf.alignment.get_alignment_length()

        # go through each seq item
        for n in range(con_len):
            # n = 1077
            # keep track of the counts of the different atoms we get
            atom_dict = {}
            num_atoms = 0

            for record in aln:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                if n < len(record.seq):
                    if record.seq[n] not in atom_dict:
                        atom_dict[record.seq[n]] = 1
                    else:
                        atom_dict[record.seq[n]] += 1

                    num_atoms += 1
            # Convert to uppercase
            atom_dict_new = dict()
            for atom in atom_dict:
                atom_dict_new[atom.upper()] = atom_dict[atom]
            atom_dict = atom_dict_new
            

            # Convert to frequencies
            for atom in atom_dict:
                atom_dict[atom] = atom_dict[atom]/num_atoms

            
            # Split ambiguities into individual DNA
            atom_dict_add = dict()
            atom_dict_remove = []
            for atom in atom_dict:
                # atom = 'K'
                if atom.replace('U','T') in ambig_to_dna.keys():
                    new_atoms = list(ambig_to_dna[atom.replace('U', 'T')])
                    for new_atom in new_atoms:
                        split_value = atom_dict[atom]/len(new_atoms)
                        if new_atom in atom_dict_add:
                            atom_dict_add[new_atom] += split_value
                        else:
                            atom_dict_add[new_atom] = split_value
                    atom_dict_remove.append(atom)
            for atom in atom_dict_remove:
                atom_dict.pop(atom)
            for atom in atom_dict_add:
                if atom in atom_dict:
                    atom_dict[atom] += atom_dict_add[atom]
                else:
                    atom_dict[atom] = atom_dict_add[atom]
            
            # Check if gaps are most frequent
            if '-' in atom_dict.keys():
                if atom_dict['-'] > 0.5:
                    consensus += '-'
                    continue
                else:
                    atom_dict.pop('-')
            
            # Check if threshold has been exceeded, if not set to 0
            thresh = 0
            max_atoms = 0
            for atom, freq in atom_dict.items():
                if freq > max_atoms:
                    max_atoms = freq
            
            if max_atoms >= threshold:
                thresh = threshold
            
            # Select atoms exceeding the threshold
            selected_atoms = []
            for atom in atom_dict:
                if atom_dict[atom] >= thresh:
                    selected_atoms.append(atom)
            
            selstring = ''.join(sorted(selected_atoms)).upper().replace('U', 'T')
            
            consensus += dna_to_ambig[selstring]
        
        if len(consensus) != con_len:
            sys.exit("Error: consensus length does not match alignment length!")
        
        return Seq(consensus)

# Main

if __name__ == "__main__":
    
    # Get options
    #args = parser.parse_args(['-t', '0.3', '-n', 'CONS00003'])
    args = parser.parse_args()
    
    
    # Read nucleotides
    # aln = AlignIO.read("/home/thomas/work/iBioGen_postdoc/MMGdatabase/phylogeny/reftree498_project/3_alignment/12_nt_aln_consex/CONS00001_COX2.fa", 'fasta')
    # aln = AlignIO.read("/tmp/source.fa", "fasta")
    try:
        aln = AlignIO.read(sys.stdin, "fasta")
    except ValueError:
        sys.exit(0)
    # Write consensus
    SeqIO.write(SeqRecord(mr_consensus(aln, args.threshold), 
                          id=args.name, description = ''), 
                sys.stdout, "fasta")
