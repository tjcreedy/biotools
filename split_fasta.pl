#!/usr/bin/env perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;

my $input;
my $outpath;

GetOptions("input=s"	=> \$input,
	   "output=s"	=> \$outpath) or die "Error getting options\n";

my $fain = Bio::SeqIO->new(-file => $input,
			   -format => "fasta");

while(my $seq = $fain->next_seq){
	my $name = $seq->id;
	$name =~ s/\s.*$//;
	$seq->id($name);
	my $faout = Bio::SeqIO->new(-file => ">$outpath/$name.fa",
				    -format => "fasta");
	$faout->write_seq($seq);
}

exit;
