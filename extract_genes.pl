#!/usr/bin/env perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;
use File::Path qw(make_path);
use File::Basename;
use Array::Utils qw(intersect);
use List::MoreUtils qw(uniq);

my @gbpaths;
my $outdir;

my $mincds = 0;
my @reqcds;

my $dontchecknames;
my $showgenes;
my $help;

my $script = basename($0, ());

my %genesort = (ATP6 => "ATP6",
		ATP8 => "ATP8",
		COX1 => "COX1",
		CO1  => "COX1",
		COI  => "COX1",
		COX2 => "COX2",
		CO2  => "COX2",
		COII => "COX2",
		COX3 => "COX3",
		CO3  => "COX3",
		COIII=> "COX3",
		CYTB => "CYTB",
		COB  => "CYTB",
		ND1  => "ND1",
		NAD1 => "ND1",
		ND2  => "ND2",
		NAD2 => "ND2",
		ND3  => "ND3",
		NAD3 => "ND3",
		ND4  => "ND4",
		NAD4 => "ND4",
		ND4L => "ND4L",
		NAD4L=> "ND4L",
		ND5  => "ND5",
		NAD5 => "ND5",
		ND6  => "ND6",
		NAD6 => "ND6");

########################################################
# USAGE
#
my $usage =<<USAGE;

Description:
	This script finds and extracts the sequences corresponding to all CDS annotations from a set of sequences from one or more genbank-format flat files. 
	
	The gene of each CDS is identified based on the /gene tag in the genbank-format flat file. Two methods are used to ensure that naming variants are correctly identified as the same gene:
		1. the gene name is converted to uppercase (so that "atp8" is the same as "ATP8")
		2. the script removes any semicolons, underscores, hyphens or spaces, and any characters following these, from gene names ("ATP8-0" -> "ATP8", "ATP8; atp8" -> "ATP8", etc)
	
	Additionally, common naming variants of mitochondrial CDS regions are hardcoded into the script (e.g. CDS regions COX1 and COI will both be output to COX1). The variants that the script knows about can be reviewed using the -showgenes option.
	
	A separate fasta will be written for each CDS gene found, containing all of the sequences of that CDS from across all the different sequences in the genbank file(s). The LOCUS name of the sequence in the genbank-format flat file will be used as the sequence header.
	
	Optionally, you can specify that CDS sequences will only be output for sequences where a minimum number of CDS annotations are met using -mincds. Similarly, you can optionally specify that certain CDS regions must be present in a sequence for it to output any CDS regions using the -reqcds option. For example, if you only wanted to retrieve sequences for CDS regions from genomes that contain COX1 and COB, you would specify -reqcds COX1 COB. The names used must conform with the standard set of output names (the second column when using the -showgenes option).
	
	By default, the script will check that sequence names fit standard naming conventions (a string of letters followed by a string of numbers, with no other characters, less than 16 characters), ignoring any input sequences that do not pass. This can be turned off.
	
Usage:

	perl $script -genbank <in.gb> [<in2.gb>] -out <out/> [-dontchecknames] [-mincds <n>] [-reqcds <CDS1> <CDS2>]
	perl $script -showgenes

Arguments:

	       genbank:  The path to one or more genbank files from which CDS regions should be extracted
	           out:  The path to a directory in which to place a fasta file for each CDS
	        mincds:  A minimum number of CDS regions that must be present in a sequence for it to output anything
	        reqcds:  A list of CDS regions that must be present in a sequence for it to output anything
	dontchecknames:  Turns off checking that sequence names fit the standard format
	     showgenes:  Prints the hardcoded conversions for CDS naming variants
	          help:  Prints out this helpful message

USAGE
#
######################################################


GetOptions("genbank=s{1,}"	=> \@gbpaths,
	   "out=s"		=> \$outdir,
	   "mincds=i"		=> \$mincds,
	   "reqcds=s{0,}"	=> \@reqcds,
	   "dontchecknames"	=> \$dontchecknames,
	   "showgenes"		=> \$showgenes,
	   "help"		=> \$help) or die "\nError getting options\n\n";

print $usage and exit if $help;
print "\t\'VARIANT\' => \'OUTPUT\'\n", Dumper \%genesort and exit if $showgenes;


# Check all values of reqcds make sense
my @knowngenes = uniq (values %genesort);
die "Error: one or more CDS names given to -reqcds don't match any known gene names\n" unless scalar (intersect(@reqcds, @knowngenes)) == scalar @reqcds;

# Make output directory
make_path($outdir);

# Set up fasta objects for outputs
my %faobjs;

# Set up hash for unrecognised gene names
my %unrec_genes;

# Work through genbank files
foreach my $gbp (@gbpaths){
	
	# Initialise genbank read object
	my $gb_in = Bio::SeqIO->new(-file => $gbp,
				    -format => "genbank",
				    -verbose => -1);
	
	# Work through sequences in object
	while(my $seq = $gb_in->next_seq){
		$seq->verbose(-1);
		my $seqname = $seq->display_id;
		#printf "Working on sequence %s from $gbp\n", $seq->display_id;
		warn "Sequence $seqname in $gbp does not fit standard naming conventions, it will be skipped\n" and next unless ($seqname =~ /^[A-Za-z_]+\d+$/ and length($seqname) <= 16 ) or $dontchecknames;
		
		# Get all CDS features from sequence
		my @cds_feats = grep {$_->primary_tag eq 'CDS' and $_->has_tag('gene')} ($seq->get_SeqFeatures);
		my $ncds_feats = scalar @cds_feats;
		
		# Check that minimum number of CDS regions are present
		warn "Sequence $seqname in $gbp has $ncds_feats CDS features, it will be skipped\n" and next unless $ncds_feats >= $mincds;
		
		# If requires any CDS regions
		if(@reqcds){
			# Extract and standardise gene names
			my @cds_feats_genes = map {$_->get_tag_values('gene')} @cds_feats;
			@cds_feats_genes = uniq @cds_feats_genes;
			@cds_feats_genes = map {$_ =~ s/[;\-_\s].*$//; uc $_} @cds_feats_genes;
			@cds_feats_genes = map {$genesort{$_} if $genesort{$_}} @cds_feats_genes;
			# Check that all required CDS are present
			warn "Sequence $seqname in $gbp does not include all required CDS regions, it will be skipped\n" and next if scalar (intersect(@reqcds, @cds_feats_genes)) < scalar @reqcds;
		}
		
		# Work through each feature
		foreach my $feat (@cds_feats){
			$feat->verbose(-1);
			my ($gene) = $feat->get_tag_values('gene');
			$gene = uc $gene;
			$gene =~ s/[;\-_\s].*$//;
			
			# Overwrite gene name if it's in the hash
			if($genesort{$gene}){
				$gene = $genesort{$gene};
			} else {
				$unrec_genes{$gene} = 1;
				warn "Found unrecognised gene name $gene in sequence $seqname in $gbp, will form a separate output\n";
			}
			
			# Check if we need to make an output object for this gene
			if(! $faobjs{$gene}){
				# Make a new object
				$faobjs{$gene} = Bio::SeqIO->new(-file => ">${outdir}/$gene.fa",
								 -format => "fasta");
			}
			
			# Compile output sequence object
			my $outseq = $feat->seq;
			$outseq->display_id($seq->display_id);
			$outseq->verbose(-1);
			# Write out sequence object
			$faobjs{$gene}->write_seq($outseq);
		}
	}
}

print "Excess files due to unrecognised gene names:\n", join ", ", keys %unrec_genes, "\n" if scalar keys %unrec_genes > 0;

exit;
