#!/usr/bin/env perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;
use List::Util qw/none/;
use List::MoreUtils qw(first_index uniq);

my %classifierpath;
my $treepath;
my $outpath;
my @taxa;
my $stringpos;
my %taxsuffix = (superfamily	=> "oidea",
		epifamily	=> "oidae",
		family		=> "[^o]idae",
		subfamily	=> "inae",
		infrafamily	=> "odd",
		tribe		=> "ini",
		subtribe	=> "ina[^e]",
		infratribe	=> "ad");

my @taxlevels = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species');

my $showconfidence;
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
	This script renames the otu terminals of newick-format trees using taxonomy data from taxonomic classifiers. Currently, it can read outputs from -rdp, -sintax and -megan. At least one of these should be supplied. Confidence levels (for RDP and SINTAX) will be shown if -showconfidence is used.
	
	The terminals of the newick tree should be named using exactly the same format as in the classifier output(s).
	
	The taxonomic levels to be used should be passed as a space-separated lower case string to -taxa, e.g. -taxa order family. For RDP and SINTAX, only the 7 major taxonomic levels will be extracted. For MEGAN, any taxonomic level between superfamily and infratribe will be extracted (using the standard latin suffixes for animal taxa). For MEGAN only, a set of consecutive taxonomic levels can additionally or alternatively be specified by passing a comma-separated set of integers defining positions along the taxonomic string to the -stringpos argument. See the example usage below
	
	If a taxonomic level cannot be extracted for a classifier, it will be ignored, silently.
	
	Multiple trees can be renamed, provided they are provided as multiple lines in a single file passed to -tree
	
	The parts of the new name will be concatenated with the character(s) specified by the -sepchar option (default ~) and surrounded by the character specified by the -bracket option (default '). Spaces can be replaced with underscores using the -replacespace option, the default is to leave spaces.
	
	Any tree terminals not in the supplied classifier output files will be silently ignored and unchanged.
	
Usage:

	perl $script -tree tree.tre -out out.tre [-sintax path -megan path -rdp path] [-taxa taxon1 taxon2 ...] [-stringpos n ] [-showconfidence -sepchar x -bracket x -replacespace]
	
	Specifying class to family on the taxonomic string \"kingdom;phylum;class;order;family;genus;species\" in MEGAN:
	perl $script -tree tree.tre -out -out.tree -megan file.txt -stringpos 3,4,5
	
	Specifying phylum, order and family in MEGAN with the same string above
	perl $script -tree tree.tre -out -out.tree -megan file.txt -stringpos 2,5 -taxa family

Arguments:

	          tree:  The path to a text file containing one or more trees
	        sintax:  The path to a SINTAX-formatted classification of sequences, one per line
	         megan:  The path to a MEGAN-formatted classification of sequences, one per line (MEGAN output readName_to_taxonPath)
	           rdp:  The path to a RDP-formatted classification of sequences, one per line
	           out:  The path to write the renamed trees
	          taxa:  A space-separated list of taxa to use for renaming
	     stringpos:  One or two integers defining the taxonomic level or set of levels to return from MEGAN
	showconfidence:  Include the confidence levels for the classifications
	       sepchar:  The character used to separate parts of the new name
	       bracket:  The character used to flank the new name
	  replacespace:  Replace spaces in new names with underscores
	          help:  Prints out this helpful message

USAGE
#
######################################################


GetOptions("tree=s"		=> \$treepath,
	   "out=s"		=> \$outpath,
	   "sintax=s"		=> \$classifierpath{SINTAX},
	   "megan=s"		=> \$classifierpath{MEGAN},
	   "rdp=s"		=> \$classifierpath{RDP},
	   "taxa=s{1,}"		=> \@taxa,
	   "stringpos=s"	=> \$stringpos,
	   "showconfidence"	=> \$showconfidence,
	   "sepchar=s"		=> \$sepchar,
	   "bracket=s"		=> \$bracket,
	   "replacespace"	=> \$replacespace,
	   "help"		=> \$help) or die "\nError getting options\n\n";

print $usage and exit if $help;

my %taxstrings;

my @classifiers = map { $_ if defined($classifierpath{$_}) } sort keys %classifierpath;

die "Error: no classification data supplied\n" if none { defined($_) } values %classifierpath;

die "Error: no tree supplied\n" unless $treepath;
die "Error: no taxa supplied\n" unless @taxa;

foreach my $classifier (@classifiers){
	open my $in_fh, '<', $classifierpath{$classifier} or die "Error opening $classifierpath{$classifier}";
	#print "-------------$classifier-----------------\n";
	while(my $line = <$in_fh>){
		chomp $line;
		$line =~ s/"//g;
		my @row = split /[\t;,:]/, $line;
		my $terminal = $row[0];
		my @taxout;
		#print "    -> $terminal\n";
		foreach my $taxlevel (@taxa){
			my $taxon;
			if($classifier eq "MEGAN"){
				if($taxsuffix{$taxlevel}){
					my $suffix = $taxsuffix{$taxlevel};
					my @taxaresults = grep(/$suffix *$/, @row);
					$taxon = $taxaresults[0] if @taxaresults;
				}
			} else {
				my $i = first_index { $_ eq $taxlevel } @taxlevels;
				if($i eq -1){
					undef $i;
				} elsif($classifier eq "SINTAX"){
					$i = ($i*2)+2;
				} elsif($classifier eq "RDP"){
					$i = ($i*3)+5;
				}
				$taxon = $row[$i] if $i;
				$taxon =~ s/\(.+\)// if $taxon and !$showconfidence;
				$taxon .= "(".$row[$i-1].")" if $i and $classifier eq "RDP" and $showconfidence;
			}
			push @taxout, $taxon if $taxon;
			
		}
		
		if($classifier eq "MEGAN" and $stringpos){
			my @stringis = split(/,/,$stringpos);
			my @taxis;
			if (@taxout){
				foreach my $taxon (@taxout){
					push @taxis, first_index{ $_ eq $taxon } @row;
				}
			}
			my @is = uniq (@stringis, @taxis);
			undef @taxout;
			foreach my $i (sort { $a <=> $b } @is){
				push @taxout, $row[$i] if defined($row[$i]);
			}
			
			#print "stringis:".join(",",@stringis)." taxis:".join(",",@taxis)." is:".join(",",@is)." taxout:".join(",",@taxout)."\n";
		}
		
		if(@taxout){
			my $outline = "$classifier:".join(',', @taxout);
			push( @{ $taxstrings{$terminal} }, $outline );
		}
	}
	
	close $in_fh;
}

my %conversion;

foreach my $terminal (keys %taxstrings){
	$conversion{$terminal} = $bracket.join($sepchar, ( ($terminal), @{ $taxstrings{$terminal} } ) ).$bracket;
	$conversion{$terminal} =~ s/ /_/g if $replacespace;
}


open my $treein, '<', $treepath or die "Error opening $treepath\n";
open my $treeout, '>', $outpath or die "Error opening $outpath for writing\n";

while(my $tree = <$treein>){
	chomp $tree;
	foreach my $seqname (keys %conversion){
		if($tree =~ /$seqname/){
			print "Renaming $seqname -> $conversion{$seqname}\n";
			$tree =~ s/(?<=[\(,])$seqname:/$conversion{$seqname}:/;
		}
	}
	print $treeout "$tree\n";
}

close $treein;
close $treeout;

exit;
