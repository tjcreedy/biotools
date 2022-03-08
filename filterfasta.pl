#!/usr/bin/env perl
#Version: 

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $fastapath;
my $filterpath;
my $taxon;
my $outpath;
my $rev;
my $reportpath;

GetOptions("fasta=s"		=> \$fastapath,
	   "filter=s"		=> \$filterpath,
	   "taxon=s"		=> \$taxon,
	   "remove"		=> \$rev,
	   "report=s"		=> \$reportpath,
	   "output=s"		=> \$outpath) or die "Error getting options\n";

# Read in filter
my %filter;

if($filterpath){
	open my $fl, '<', $filterpath or die "Error reading $filterpath\n";
	my $frow = <$fl>;
	chomp $frow;
	if($frow =~ /^\(/){
		%filter = map { $_ => 1 } ($frow =~ /[\(,]([^:\(]+):/g);
	} else {
		seek $fl, 0, 0;
		while(my $row = <$fl>){
			chomp $row;
			$row =~ s/"//g;
			my @line = split /[\t,]/,$row;
			
			if($taxon and $line[1] =~ /$taxon/){
				$filter{$line[0]} = 1;
			} elsif(!$taxon){
				$filter{$line[0]} = 1;
			}
		}
	}
	
	close $fl;
}

#print Dumper \%filter;

#Read and write fasta
open my $fo, '>', $outpath or die "Error opening $outpath for writing\n";
open my $fi, '<', $fastapath or die "Error reading $fastapath\n";
my $printswitch;
while(my $row = <$fi>){
	chomp $row;
	if($row =~ /^>(.+)$/){
		if($filterpath){
			if($filter{$1}){
				$filter{$1}++;
				$printswitch = $rev ? 0 : 1;
			} else {
				$printswitch = $rev ? 1 : 0;
			}
		} elsif($taxon){
			my $head = $1;
			if($head =~ /$taxon/){
				$filter{$taxon}++;
				$printswitch = $rev ? 0 : 1;
			} else {
				$printswitch = $rev ? 1 : 0;
			}
		}
	}
	if($printswitch){
		print $fo "$row\n";
	}
}
close $fo;
close $fi;

#print Dumper \%filter;

# Report on filtering
if($reportpath){
	my %report;
	foreach my $head (keys %filter){
		if($filter{$head} < 2){
			push @{ $report{nomatch} }, $head;
		} elsif($filter{$head} == 2){
			push @{ $report{onematch} }, $head;
		} else {
			push @{ $report{twoplusmatch} }, $head;
		}
	}
	#print Dumper \%report;
	open my $ro, '>', $reportpath or die "Error opening $reportpath for writing\n";
	
	if($report{twoplusmatch}){
		print $ro "\n===== These headers were found two or more times in the fasta =====\n";
		foreach my $head (values @{ $report{twoplusmatch} }){
			print $ro "$head\n";
		}
	}
	if ($report{nomatch}){
		print $ro "\n===== These headers were not found in the fasta ======\n";
		foreach my $head (values @{ $report{nomatch} }){
			print $ro "$head\n";
		}
	}
	if ($report{onematch}){
		print $ro "\n===== These headers were found exactly once in the fasta =====\n";
		foreach my $head (values @{ $report{onematch} }){
			print $ro "$head\n";
		}
	}
	
	close $ro;
}

exit;
