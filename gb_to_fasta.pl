#!/usr/bin/perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;

my $gbpath;
my $fapath;

GetOptions("genbank=s"	=> \$gbpath,
	   "fasta=s"	=> \$fapath) or die "Error getting options\n";

my $gb = Bio::SeqIO->new(-file => $gbpath, -format => "genbank");
my $fa = Bio::SeqIO->new(-file => ">$fapath", -format => "fasta");

while(my $seq = $gb->next_seq){
	$fa->write_seq($seq);
}

exit;
