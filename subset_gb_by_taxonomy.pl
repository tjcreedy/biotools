#!/usr/bin/env perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;
use Bio::TreeIO;
use List::Util qw/shuffle/;

my $gbpath;
my $taxlevel;
my $outname;
my @outfiles = ("subset", "excluded", "no_taxon", "excess_taxon");
my $nout = 1;
my $incsmall;
my %taxsuffix = (superfamily	=> "oidea",
		epifamily	=> "oidae",
		family		=> "[^o]idae",
		subfamily	=> "inae",
		infrafamily	=> "odd",
		tribe		=> "ini",
		subtribe	=> "ina[^e]",
		infratribe	=> "ad");
my $suffix;
my $random;
my $help;

my $script = basename($0, ());

########################################################
# USAGE
#
my $usage =<<USAGE;

Description:

	This script subsets a -genbank file according to groups of the specified taxonomic level (-taxlevel), outputting the specified -number of entries for each taxonomic group (default 1). The genbank file must have taxonomy annotations otherwise the script will not work. It does not look up taxonomy from NCBI taxids. 
	
	Taxonomic levels are identified based on the standard animal latin naming convention suffixes for levels from superfamily to infratribe inclusive. Alternatively, a custom suffix can be supplied using the -customsuffix. Anything supplied to -customsuffix will be ignored if a taxonomic level is supplied. Note that suffixes must be the end of the level, and that suffixes that are subsets of any other possible suffix should use regular expressions to avoid duplicate matches, which will result in exclusions. For example, the family suffix \"idae\" is a subset of \"oidae\", so the former will match twice where the taxonomy includes an epifamily. To avoid this, the suffix [^o]idae is used internally. Custom suffixes are not checked for possible collision with standard suffixes.
	
	The script will pick the first entry or entries for each taxonomic group by default, but can optionally select randomly from all entries available using the -random option.
	
	If fewer than the requested number of entries are available for a taxonomic group, the script will either exclude all entries in the group (the default), or output all available entries using the -includesmall option.
	
	Up to four genbak-format files will be output with the prefix -outname. The _subset file will contain the selected subsets and the _excluded file will contain all other entries with an identified taxonomic group. If for any entries no taxonomic group could be identified, or multiple taxonomic groups were identified, the relevant sequences will be output to the _notaxon and _excess_taxon files respectively.

Usage:

	perl $script -genbank file.gb -outname prefix -taxlevel level [-number n -random -includesmall]
	perl $script -genbank file.gb -outname prefix -customsuffix suffix [-number n -random -includesmall]
	
	For example, to select 2 random entries for each family, including any singletons:
	perl $script -genbank file.gb -outname out -taxlevel family -number 2 -random -includesmall

Arguments:

	       genbank:  The path to a genbank-format files from which entries will be selected
	       outname:  The suffix of the output file name
	      taxlevel:  The taxonomy level to group entries by (lowercase, superfamily-infratribe)
	        number:  The number of entries to select from each taxonomic group (Default: 1)
	        random:  Select entries randomly
	  includesmall:  Include entries from a taxonomic group where the number of entries in that group is < number
	          help:  Prints out this helpful message

USAGE
#
######################################################


GetOptions("genbank=s"		=> \$gbpath,
	   "number=i"		=> \$nout,
	   "taxlevel=s"		=> \$taxlevel,
	   "customsuffix=s"	=> \$suffix,
	   "outname=s"		=> \$outname,
	   "random"		=> \$random,
	   "includesmall"	=> \$incsmall,
	   "help"		=> \$help) or die "\nError getting options\n\n";

print $usage and exit if $help;

die "\nError: please supply a genbank-format file\n" unless $gbpath;
die "\nError: please specify an output name\n" unless $outname;


if($suffix){
	$taxlevel = "custom"
} elsif($taxlevel){
	my $posslevels = join "\n\t", keys %taxsuffix;
	die "\nError, taxonomic level $taxlevel not recognised. It must be one of the following: \n\t$posslevels\n" unless $taxsuffix{$taxlevel};
	
	$suffix = $taxsuffix{$taxlevel};
} else {
	die "\nError, please enter a taxonomic level or a custom suffix\n";
}

(my $nicesuffix = $suffix) =~ s/\[.*\]//;

print "\nUsing $taxlevel suffix \"$nicesuffix\" to identify the taxonomic level to group by and subset from\n";


my $gb_in = Bio::SeqIO->new(-file => $gbpath,
			    -format => "genbank");

my %outhandle;

foreach my $of (@outfiles){
	$outhandle{$of} = Bio::SeqIO->new(-file => ">${outname}_$of.gb",
					  -format => "genbank");
}

my %matchrecords;

my $nseq;
my %outseqs;

while(my $seq = $gb_in->next_seq){
	$nseq++;
	$seq->verbose(-1);
	my @taxonomy = $seq->species->classification();
	
	my @taxa = grep(/$suffix *$/, @taxonomy);
	
	my $result;
	
	if(scalar @taxa == 1){
		my $taxon = $taxa[0];
		#if($matchrecords{$taxon}){
			push(@{ $matchrecords{$taxon} }, $seq);
		#} else {
		#	$matchrecords{$taxon} = [ $seq ];
		#}
	} elsif(scalar @taxa == 0){
		$result = "no_taxon";
	} else {
		$result = "excess_taxon";
	}
	
	if($result){
		$outhandle{$result}->write_seq($seq);
		$outseqs{$result}++;
	}
}

print "\nLoaded $nseq sequences from $gbpath\n\n";

foreach my $taxon (keys %matchrecords){
	my @seqs;
	
	if($random){
		@seqs = shuffle @{ $matchrecords{$taxon} };
	} else {
		@seqs = @{ $matchrecords{$taxon} };
	}
	delete $matchrecords{$taxon};
	
	my $output_group = (scalar @seqs >= $nout or $incsmall);
	
	my $written = 0;
	foreach my $i (0..$#seqs){
		my $result;
		if ($i < $nout and $output_group){
			$result = "subset";
			$written++;
		} else {
			$result = "excluded";
		}
		my $seq = $seqs[$i];
		$seq->verbose(-1);
		$outhandle{$result}->write_seq($seqs[$i]);
		$outseqs{$result}++;
	}
	
	printf "Written $written entries for $taxon (%d total)\n", scalar @seqs;
}

$outseqs{excluded} = 0 unless $outseqs{excluded};

if($outseqs{subset}){
	print "\nCompleted writing a subset of $outseqs{subset} entries, excluding $outseqs{excluded} entries.\n";
} else {
	print "\nNo entries were selected with the given parameters. $outseqs{excluded} entries had an identifiable taxon but were excluded\n";
}

if($outseqs{no_taxon}){
	print "\n$outseqs{no_taxon} entries did not have any taxa matching the suffix \"$nicesuffix\"\n";
} else {
	unlink "${outname}_no_taxon.gb";
}

if($outseqs{excess_taxon}){
	print "\n$outseqs{excess_taxon} entries had more than one taxa matching the suffix \"$nicesuffix\"\n";
} else {
	unlink "${outname}_excess_taxon.gb";
}

exit;
