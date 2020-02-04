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

my $gbin = Bio::SeqIO->new(-file => $input,
			   -format => "genbank");

while(my $seq = $gbin->next_seq){
	my $name = $seq->display_id;
	$name =~ s/\s.*$//;
	my $gbout = Bio::SeqIO->new(-file => ">$outpath/$name.gb",
				    -format => "genbank");
	$gbout->write_seq($seq);
}

exit;
