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

my @gbpaths;
my $outdir;

my $minregion = 1;
my @reqregion;
my @regiontypes = ('CDS');
my $keepframe;
my $organism;

my $dontchecknames;
my $showgenes;
my $help;

my $script = basename($0, ());

my $genenames =<<GENES;
ATP6;CDS:ATP SYNTHASE F0 SUBUNIT 6,ATP6,APT6,ATP SYNTHASE A0 SUBUNIT 6,ATP SYNTHASE SUBUNIT 6
ATP8;CDS:ATP SYNTHASE F0 SUBUNIT 8,ATP8,APT8,ATP SYNTHASE A0 SUBUNIT 8,ATP SYNTHASE SUBUNIT 6
COX1;CDS:CYTOCHROME C OXIDASE SUBUNIT 1,CYTOCHROME OXIDASE SUBUNIT I,CYTOCHROME C OXIDASE SUBUNIT I,COX1,COXI,CO1,COI
COX2;CDS:CYTOCHROME C OXIDASE SUBUNIT 2,CYTOCHROME OXIDASE SUBUNIT II,CYTOCHROME C OXIDASE SUBUNIT II,COX2,COXII,CO2,COII
COX3;CDS:CYTOCHROME C OXIDASE SUBUNIT 3,CYTOCHROME OXIDASE SUBUNIT III,CYTOCHROME C OXIDASE SUBUNIT III,COX3,COXII,CO3,COIII
CYTB;CDS:CYTOCHROME B,CYTB,CYB,COB
ND1;CDS:NADH DEHYDROGENASE SUBUNIT 1,ND1,NAD1,NSD1,NADH1
ND2;CDS:NADH DEHYDROGENASE SUBUNIT 2,ND2,NAD2,NSD2,NADH2
ND3;CDS:NADH DEHYDROGENASE SUBUNIT 3,ND3,NAD3,NSD3,NADH3
ND4;CDS:NADH DEHYDROGENASE SUBUNIT 4,ND4,NAD4,NSD4,NADH4
ND4L;CDS:NADH DEHYDROGENASE SUBUNIT 4L,ND4L,NAD4L,NSD4L,NADH4L
ND5;CDS:NADH DEHYDROGENASE SUBUNIT 5,ND5,NAD5,NSD5,NADH5
ND6;CDS:NADH DEHYDROGENASE SUBUNIT 6,ND6,NAD6,NSD6,NADH6
18S;rRNA:18S RIBOSOMAL RNA,18SRRN,18S
RRNL;rRNA:LARGE SUBUNIT RIBOSOMAL RNA,RRNL,L,LSU,LARGE,L-RRNA,16S RIBOSOMAL RNA,16SRRN,16S,16S-RNA
RRNS;rRNA:SMALL SUBUNIT RIBOSOMAL RNA,RRNS,S,SSU,SMALL,S-RRNA,12S RIBOSOMAL RNA,12SRRN,12S,12S-RNA
TRNA;tRNA:TRNA,TRNA (AGC),TRNA-AGC,TRNA (UGC),TRNA (CGC),TRNA-CGC,TRNA (TGC),TRNA-UGC,TRNA-TGC,TRNA-ALA
TRNC;tRNA:TRNC,TRNC (GCA),TRNC-GCA,TRNC (ACA),TRNC-ACA,TRNA-CYS
TRND;tRNA:TRND,TRND (GUC),TRND (GTC),TRND-GUC,TRND-GTC,TRND (AUC),TRND (ATC),TRND-AUC,TRND-ATC,TRNA-ASP,TRND(---)
TRNE;tRNA:TRNE,TRNE (UUC),TRNE (TTC),TRNE-UUC,TRNE-TTC,TRNE (CUC),TRNE (CTC),TRNE-CUC,TRNE-CTC,TRNA-GLU
TRNF;tRNA:TRNF,TRNF (GAA),TRNF-GAA,TRNF (AAA),TRNF-AAA,TRNA-PHE,TRNF(---)
TRNG;tRNA:TRNG,TRNG (GCC),TRNG-GCC,TRNG (UCC),TRNG (TCC),TRNG-UCC,TRNG-TCC,TRNA-GLY
TRNH;tRNA:TRNH,TRNH (GUG),TRNH (GTG),TRNH-GUG,TRNH-GTG,TRNH (AUG),TRNH (ATG),TRNH-AUG,TRNH-ATG,TRNA-HIS
TRNI;tRNA:TRNI,TRNI (GAU),TRNI (GAT),TRNI-GAU,TRNI-GAT,TRNI (AAU),TRNI (AAT),TRNI-AAU,TRNI-AAT,TRNA-ILE,TRNI(---)
TRNK;tRNA:TRNK,TRNK (CUU),TRNK (CTT),TRNK-CUU,TRNK-CTT,TRNK (UUU),TRNK (TTT),TRNK-UUU,TRNK-TTT,TRNA-LYS
TRNL;tRNA:TRNL,TRNL (UAA),TRNL (TAA),TRNL-UAA,TRNL-TAA,TRNL (UAG),TRNL (TAG),TRNL-UAG,TRNL-TAG,TRNA-LEU,TRNA-LEU2,TRNL2 (TAA),TRNL(---),TRNA-LEU1,TRNL1(CAG),TRNL2(UUR),TRNL1(AAG),TRNL1,TRNL1(TAG),TRNL(CUN),TRNL1 (TAG),TRNL2,TRNL1(CUN),TRNL2(TAA)
TRNM;tRNA:TRNM,TRNM (CAU),TRNM (CAT),TRNM-CAU,TRNM-CAT,TRNM (UAU),TRNM (TAT),TRNM-UAU,TRNM-TAT,TRNA-MET,TRNM(---)
TRNN;tRNA:TRNN,TRNN (GUU),TRNN (GTT),TRNN-GUU,TRNN-GTT,TRNA-ASN,TRNN(---)
TRNP;tRNA:TRNP,TRNP (UGG),TRNP (TGG),TRNP-UGG,TRNP-TGG,TRNA-PRO
TRNQ;tRNA:TRNQ,TRNQ (UUG),TRNQ (TTG),TRNQ-UUG,TRNQ-TTG,TRNQ (CUG),TRNQ (CTG),TRNQ-CUG,TRNQ-CTG,TRNA-GLN,TRNQ(A--),TRNQ(---)
TRNR;tRNA:TRNR,TRNR (CCG),TRNR-CCG,TRNR (GCG),TRNR-GCG,TRNR (UCG),TRNR (TCG),TRNR-UCG,TRNR-TCG,TRNA-ARG
TRNS;tRNA:TRNS,TRNS (AGA),TRNS-AGA,TRNS (UCU),TRNS (TCT),TRNS-UCU,TRNS-TCT,TRNS (UGA),TRNS (TGA),TRNS-UGA,TRNS-TGA,TRNS (GCU),TRNS (GTU),TRNS-GCU,TRNL-GTU,TRNA-SER,TRNA-SER1,TRNA-SER2,TRNS2 (TGA),TRNS2(UCN),TRNS1,TRNS1 (TCT),TRNS1(GCT),TRNS1(AGN),TRNS2(TGA),TRNS1(ACT),TRNS(---),TRNS2,TRNS2(CGA),TRNS1(TCT)
TRNT;tRNA:TRNT,TRNT (UGU),TRNT (TGT),TRNT-UGU,TRNT-TGT,TRNA-THR,TRNT(---)
TRNV;tRNA:TRNV,TRNV (GAC),TRNV-GAC,TRTV (CAC),TRTV-CAC,TRNV (UAC),TRNV (TAC),TRNV-UAC,TRNV-TAC,TRNA-VAL,TRNV(CAC),TRNV(AAC)
TRNW;tRNA:TRNW,TRNW (UCA),TRNW (TCA),TRNW-UCA,TRNW-TCA,TRNA-TRP,TRNW(-CA)
TRNY;tRNA:TRNY,TRNY (GUA),TRNY (GTA),TRNY-GUA,TRNY-GTA,TRNY (AUA),TRNY (ATA),TRNY-AUA,TRNY-ATA,TRNY (UCA),TRNY (TCA),TRNY-UCA,TRNY-TCA,TRNA-TYR,TRNY(---)
TRNX;tRNA:TRNX,TRNX (CUA),TRNX (CTA),TRNX-CUA,TRNX-CTA
GENES

my ($known_genes, $genesort, $generegion) = sort_genes($genenames);

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
		my $species = $seq->species->node_name;
		$species =~ tr/ /_/;
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
				
				# Add to unrecognised genes if unrecognised
				push @{$unrec_genes{$gene}}, $seqname unless ${$genesort}{$gene};
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
	print "Excess files due to unrecognised gene names:\n", join (", ", keys %unrec_genes), "\n\nSource sequences of unrecognised gene names:\n" ;
	
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
		my @values = split("[;:,]", $line);
		push @genenames, $values[0];
		$generegion{$values[0]} = $values[1];
		foreach my $var (@values[2..$#values]){
			$genesort{$var} = $values[0];
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
	$name =~ s/\//_/;
	return($name)
}

