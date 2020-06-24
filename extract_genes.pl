#!/usr/bin/env perl
#By Thomas J Creedy, thomas@tjcreedy.co.uk 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;
use File::Path qw(make_path);
use File::Basename;
use Array::Utils qw(intersect unique);
use List::MoreUtils qw(uniq);
use LWP::Simple;

my @gbpaths;
my $outdir;

my $minregion = 1;
my @reqregion;
my @regiontypes = ('CDS');
my $keepframe;
my $organism;
my $writeunknowns;

my $dontchecknames;
my $showgenes;
my $help;

my $script = basename($0, ());

my $genenames = get 'https://raw.githubusercontent.com/tjcreedy/biotools/master/gene_name_variants.txt';

my ($known_genes, $genesort, $generegion) = sort_genes($genenames);

print Dumper $genesort;
exit;

########################################################
# USAGE
#
my $usage =<<USAGE;

Description:
	This script finds and extracts sequences corresponding to CDS, rRNA and/or tRNA annotations from a set of sequences from one or more genbank-format flat files. Specify which regions to extract using the -regiontypes argument (default is -regiontypes CDS rRNA)
	
	The region name is identified based on the /gene, /label or /product tag in the genbank-format flat file. Two methods are used to ensure that naming variants are correctly identified as the same gene:
		1. the gene name is converted to uppercase (so that "atp8" is the same as "ATP8")
		2. the script removes any semicolons, underscores, hyphens or spaces, and any characters following these, from gene names ("ATP8-0" -> "ATP8", "ATP8; atp8" -> "ATP8", etc) - except if the first part of the name is "MT-".
	
	Additionally, common naming variants of mitochondrial regions are hardcoded into the script (e.g. CDS regions COX1 and COI will both be output to COX1). The variants that the script knows about can be reviewed using the -showgenes option.
	
	A separate fasta will be written for each region found, containing all of the sequences of that CDS from across all the different sequences in the genbank file(s). The LOCUS name of the sequence in the genbank-format flat file will be used as the sequence header. Alternatively, you can try to use use the organism field of the genbank file as the sequence header using --organism: note that this is very likely to give blank or non-unique headers, and if the script fails to parse the species name it will revert to the LOCUS field.
	
	Optionally, you can specify that sequences will only be output for entries where a minimum number of regions are met using -minregion. Similarly, you can optionally specify that certain regions must be present in a sequence for it to output any regions using the -reqregion option. For example, if you only wanted to retrieve sequences for regions from genomes that contain COX1 and COB, you would specify -reqregion COX1 COB. The names used must conform with the standard set of output names (the second column when using the -showgenes option). All available regions will still be output.
	
	In rare cases, an annotation is truncated by the end of a contig at the 5' end, which may leave an incomplete codon at the start of an output sequence and cause downstream issues. To detect this, the script can optionally remove these excess 1 or 2 bases using the -keepframe option. Note this only works if the /codon_start tag is present in the genbank-format features. 

Usage:

	perl $script -genbank <in.gb> [<in2.gb>] -out <out/> [-minregion <n>] [-reqregion <CDS1> <CDS2>] [-regiontypes CDS tRNA rRNA] [-keepframe]
	perl $script -showgenes

Arguments:

	       genbank:  The path to one or more genbank files from which CDS regions should be extracted
	           out:  The path to a directory in which to place a fasta file for each CDS
	     minregion:  A minimum number of regions that must be present in a sequence for it to output anything
	     reqregion:  A list of regions that must be present in a sequence for it to output anything
	   regiontypes:  One or more region types to extract, currently CDS rRNA and/or tRNA
	     showgenes:  Prints the hardcoded conversions for region naming variants
	      organism:  Prints the organism name, if retreivable, instead of the LOCUS
	 writeunknowns:  Writes out fastas for unidentifiable annotations anyway
	     keepframe:  Removes excess out-of-frame bases at the beginning of truncated annotations
	          help:  Prints out this helpful message

USAGE
#
######################################################


GetOptions("genbank=s{1,}"	=> \@gbpaths,
	   "out=s"		=> \$outdir,
	   "minregion=i"	=> \$minregion,
	   "reqregion=s{0,}"	=> \@reqregion,
	   "regiontypes=s{1,3}"	=> \@regiontypes,
	   "showgenes"		=> \$showgenes,
	   "organism"		=> \$organism,
	   "writeunknowns"	=> \$writeunknowns,
	   "keepframe"		=> \$keepframe,
	   "help"		=> \$help) or die "\nError getting options\n\n";

print $usage and exit if $help;

print "\t\'VARIANT\' => \'OUTPUT\'\n", Dumper $genesort and exit if $showgenes;

die "Error: please supply the path(s) to one or more genbank files and the path to an output directory\n" unless @gbpaths and $outdir;

# Check all values of reqregion make sense
die "Error: one or more region names given to -reqregion don't match any known gene names\n" unless scalar (intersect(@reqregion, @$known_genes)) == scalar @reqregion;

# Check all values of regiontypes make sense
my @alltypes = ('CDS', 'rRNA', 'tRNA');
die "Error: supply one or more region types to -regiontypes\n" unless scalar @regiontypes > 0;
die "Error: one or more region types given to -regiontypes don't match any known region types\n" unless scalar (unique(@regiontypes, @alltypes)) == 3;

my %regtypes = map {$_ => 1} @regiontypes;

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
		
		# Extract sequence name and species name
		my $seqname = $seq->display_id;
		my $species;
		if($organism){
			$species = $seq->species->node_name;
			$species =~ tr/ /_/;
		}
		# Get all features with sufficient identification from object
		my @all_feats = grep {$_->has_tag('gene') or $_->has_tag('product') or $_->has_tag('label')} ($seq->get_SeqFeatures);
		
		# Set up containers
		my %found_genes;
		my $checked_features = 0;
		
		#printf "Seqname = $seqname, Feats = %d, max_checked_feats = %d\n", $#all_feats, $#all_feats * 2 + 1;
		
		# Work through features to extract data
		while($checked_features < ($#all_feats * 2 + 1) and scalar keys %found_genes < scalar @all_feats){
			
			# Assess current position in search
			my $round = $checked_features <= $#all_feats ? 0 : 1;
			my $i = $round ? $checked_features - $#all_feats - 1 : $checked_features;
			
			# Extract feature
			my $feat = $all_feats[$i];
			$feat->verbose(-1);
			
			my $feattype = $feat->primary_tag;
			
			# Get clean gene name and convert
			my $gene = get_clean_name($feat);
			$gene = ${$genesort}{$gene} if ${$genesort}{$gene};
			
			#print "Checked_features = $checked_features, Round = $round, i = $i, feattype = $feattype, gene = $gene ...";
			
			
			# If gene has not already been extracted, it is the first round of checks and this is a required type of feature or it is the second round of checks and this is a gene feature and the gene type matches those required
			
			if( not $found_genes{$gene} and 
					( $regtypes{$feattype} or (
						$feat->primary_tag eq 'gene' and 
						$round and 
						${$generegion}{$gene} and 
						$regtypes{${$generegion}{$gene}}
						)
					)){
				#print "extracting\n";
				
				# If unrecognised
				if( not ${$genesort}{$gene}){
					push @{$unrec_genes{$gene}}, $seqname;
					$checked_features++ and next unless $writeunknowns;
				}
				
				# Make a new output object for this gene if needed
				if(! $faobjs{$gene}){
					$faobjs{$gene} = Bio::SeqIO->new(-file => ">${outdir}/$gene.fa",
									 -format => "fasta");
				}
				
				# Check if gene sequence is truncated
				my $codon_start = 1;
				if($feat->has_tag('codon_start') and $keepframe){
					my @values = $feat->get_tag_values("codon_start");
					die "Error: $seqname $feattype annotation $gene has multiple \\codon_start values!\n" if(scalar @values > 1);
					$codon_start = $values[0]
				}
				
				# Compile output sequence object and store
				my $outseq = $feat->seq;
				$outseq = $outseq->trunc($codon_start, $outseq->length());
				if($organism){
					$outseq->display_id($species) or $outseq->display_id($seqname);
				} else {
					$outseq->display_id($seqname);
				}
				$outseq->description("");
				$outseq->verbose(-1);
				
				$found_genes{$gene} = $outseq;
				
			} else {
				#print "skipping\n";
			}
			$checked_features++;
		}
		
		# Run checks against threshold number and content of found genes
		my @found_gene_names = keys %found_genes;
		my $n_found_known_genes = scalar (intersect( @found_gene_names, @{$known_genes}));
		warn "Sequence $seqname in $gbp has $n_found_known_genes known features with sufficient information, it will be skipped\n" and next unless $n_found_known_genes >= $minregion;
		
		warn "Sequence $seqname in $gbp does not include all required regions, it will be skipped\n" and next if scalar (intersect(@reqregion, @found_gene_names)) < scalar @reqregion and @reqregion;
		
		# Write out found genes
		$faobjs{$_}->write_seq($found_genes{$_}) foreach(keys %found_genes);
		
	}
}

if(scalar keys %unrec_genes > 0){
	print "Unrecognised gene names:\n", join (", ", keys %unrec_genes), "\n\nSource sequences of unrecognised gene names:\n" ;
	
	foreach my $ugene (keys %unrec_genes){
		print "$ugene: ", join (", ", @{$unrec_genes{$ugene}}), "\n";
	}
}

exit;

sub sort_genes{
	my ($genenametxt) = @_;
	my @genenames;
	my %genesort;
	my %generegion;
	foreach my $line (split("\n", $genenametxt)){
		my ($name, $region, $product, @values) = split("[;:,]", $line);
		push @values, ($name, $product);
		print Dumper \@values;
		push @genenames, $name;
		$generegion{$name} = $product;
		foreach my $var (@values){
			$genesort{$var} = $name;
		}
	}
	return(\@genenames, \%genesort, \%generegion)
}

sub get_clean_name{
	my ($feat) = @_;
	my $name;
	
	if($feat->has_tag('gene')){
		($name) = $feat->get_tag_values('gene');
	} elsif($feat->has_tag('label')){
		($name) = $feat->get_tag_values('label');
	}
	if($name){
		$name = uc $name;
		$name =~ s/[(;_ ].*$//;
		$name =~ s/(?<!MT)-.*$//;
	} elsif($feat->has_tag('product')){
		($name) = $feat->get_tag_values('product');
		$name = uc $name;
	} else {
		die "Error: feature passed to get_clean_name subroutine without gene, label or product tag\n";
	}
	$name =~ s/[\/,\\:;]/_/g;
	return($name)
}

