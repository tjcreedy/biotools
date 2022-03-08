# biotools
This repository comprises various standalone scripts that are useful for manipulating DNA sequence 
data. These are summarised here, see the documentation of each script using the help argument.

**annotation_distances.py** reports the number of bases between each consecutive pair of 
annotations for each record in a genbank format file

**translate.py** translates nucleotide sequences in a fasta into amino acids

**backtranslate.py** takes as input a aligned amino acid sequences and their corresponding
unaligned nucleotide sequences and returns a nucleotide alignment

**consensus.py** outputs the majority rule consensus of an input alignment

**extract_genes.pl** extracts the nucleotide sequences corresponding to selected annotations from 
the records of a genbank file

**filter_fasta_by_fasta.py** filters the sequences in a fasta based on their presence in a second 
fasta file

**filter_fasta_by_sintax.py** filters the sequences in a fasta according to the taxonomy in 
a sintax-format file

**filterfasta.pl** filters the sequences in a fasta according to a text file of required headers 
or a specified search term

**fixgb.py** corrects genbank file LOCUS lines to work properly in biopython

**gb_to_fasta.pl** extracts the sequences from a genbank file and outputs them as a fasta file

**get_genbanks.py** downloads genbank-format files from NCBI GenBank by accession number

**get_NCBI_Taxonomy.py** retrieves taxonomy information from the NCBI Taxonomy database by taxid 
numbers


**MBC_refmatcher.py** integrates the results of searching an OTU/ASV file against a local set of 
references into the OTU/ASV file and a read mapping table

**muldemux.py** wraps cutadapt over many files and indices

**reducealign.py** removes redundant columns from an alignment

**split_fasta.pl** splits a multifasta into individual fastas for each sequence

**split_fasta_by_label.py** splits a multifasta into separate fastas according to labels in the 
sequence headers

**split_genbank.pl** splits a genbank-format file with multiple records into individual files for 
each record

**subset_fasta.py** randomly subsets a given proportion or number of sequences from a fasta

**subset_gb_by_taxonomy.pl** extracts sequences from a genbank-format file according to taxonomy

**trimalign.py** removes sequences from an alignment according to regular expressions, optionally 
dropping gap-only columns afterwards.

## Deprecated scripts

**isPCR.py** performs a very basic attempt at in silico PCR - cutadapt is better than this script

**rename_newick_with_classifiers.pl** renames a phylogeny using a table - phylabel.R in the 
(phylostuff)[https://github.com/tjcreedy/phylostuff] repository is better than this
**rename_newick_with_gb.pl** renames a phylogeny using metadata parsed from a genbank file



