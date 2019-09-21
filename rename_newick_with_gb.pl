#!/usr/bin/env perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;
use Bio::TreeIO;
use File::Basename;

my @gbpaths;
my $treepath;
my $outpath;

my %fieldmethods = ('LOCUS'		=> 'display_id',
		    'DEFINITION'	=> 'description',
		    'ACCESSION'		=> 'accession_number',
		    'VERSION'		=> 'version',
		    'LENGTH'		=> 'length');
my @fieldstrings;
my $taxitems;
my $replacespace;
my $sepchar = "~";
my $bracket = "\'";
my $help;

my $script = basename($0, ());

########################################################
# USAGE
#
my $usage =<<USAGE;

Description:
	This script renames the terminals of newick-format trees using data from a GenBank-format file by matching the current names of the tree against the LOCUS field of the genbank file.
	
	Multiple GenBank-format files can be passed to the -genbank argument. Any entries in the GenBank file not found in the tree(s) are silently ignored.
	
	Multiple trees can be renamed, provided they are provided as multiple lines in a single file passed to -tree
	
	Specify one or more of the following fields in a space separated list to the -fields argument:
		LOCUS
		DEFINITION
		ACCESSION
		VERSION
		LENGTH
		ORGANISM
		TAXONOMY
		NCDS
	
	The TAXONOMY field corresponds to the semicolon-separated string given immediately after the ORGANISM field in a GenBank-format file, plus the ORGANISM. Semicolons will be replaced with commas. The number of taxonomic levels to be returned can be controlled using the -ntax option, starting at the lowest taxonomic level. Thus the use of -f ORGANISM and -f TAXONOMY -ntax 1 is exactly identical.
	
	The NCDS field counts the number of CDS features present in the GenBank entry.
	
	The selected fields will be concatenated with the character(s) specified by the -sepchar option (default ~) and surrounded by the character specified by the -bracket option (default '). Spaces can be replaced with underscores using the -replacespace option, the default is to leave spaces.
	
	Any tree terminals not in the supplied GenBank-format files will be silently ignored and unchanged.
	
Usage:

	perl $script -genbank <in.gb> [<in2.gb>] -tree <tree.tre> -out <out.tre> -fieldstring FIELD [FIELD] [-taxitems n] [-replacespace]

Arguments:

	       genbank:  The path to one or more genbank files from which metadata will be extracted
	          tree:  The path to a text file containing one or more trees
	           out:  The path to write the renamed trees
	        fields:  A space-separated list of fields to use for renaming
	          ntax:  If TAXONOMY is specified, the number of taxonomic levels to include
	  replacespace:  Replace spaces in new names with underscores
	          help:  Prints out this helpful message

USAGE
#
######################################################


GetOptions("genbank=s{1,}"	=> \@gbpaths,
	   "tree=s"		=> \$treepath,
	   "out=s"		=> \$outpath,
	   "fieldstring=s{1,}"	=> \@fieldstrings,
	   "ntax=i"		=> \$taxitems,
	   "sepchar=s"		=> \$sepchar,
	   "bracket=s"		=> \$bracket,
	   "replacespace"	=> \$replacespace,
	   "help"		=> \$help) or die "\nError getting options\n\n";

print $usage and exit if $help;

my %conversion;

foreach my $gbpath (@gbpaths){
	my $gb_in = Bio::SeqIO->new(-file => $gbpath,
				    -format => "genbank");
	
	while(my $seq = $gb_in->next_seq){
		
		my $seqname = $seq->display_id();
		my @fieldarray;
		
		die "Error: $seqname in $gbpath is a duplicate of a previous entry!\n" if $conversion{$seqname};
		
		foreach my $f (@fieldstrings){
			my $outitem;
			if($fieldmethods{$f}){
				my $method = $fieldmethods{$f};
				$outitem = $seq->$method;
			} elsif($f eq "ORGANISM"){
				$outitem = $seq->species->node_name();
			} elsif($f eq "TAXONOMY"){
				my @taxa = $seq->species->classification();
				$taxitems = scalar @taxa if(!$taxitems or $taxitems > scalar @taxa);
				$outitem =  join(";", reverse(@taxa[0..$taxitems-1]));
				$outitem =~ s/;/,/g;
			} elsif($f eq "NCDS"){
				my @cds_feats = grep {$_->primary_tag eq 'CDS' and $_->has_tag('gene')} ($seq->get_SeqFeatures);
				$outitem = scalar @cds_feats;
			}
			$outitem =~ s/\.$//;
			push @fieldarray, $outitem;
		}
		
		my $outline = join($sepchar, @fieldarray);
		$outline =~ s/ /_/g if $replacespace;
		$outline = "$bracket$outline$bracket";
		$conversion{$seqname} = $outline;
	}
}

open my $treein, '<', $treepath or die "Error opening $treepath\n";
open my $treeout, '>', $outpath or die "Error opening $outpath for writing\n";

while(my $tree = <$treein>){
	chomp $tree;
	foreach my $seqname (keys %conversion){
		if($tree =~ /$seqname/){
			print "Renaming $seqname -> $conversion{$seqname}\n";
			$tree =~ s/$seqname/$conversion{$seqname}/;
		}
	}
	print $treeout "$tree\n";
}

close $treein;
close $treeout;

exit;
