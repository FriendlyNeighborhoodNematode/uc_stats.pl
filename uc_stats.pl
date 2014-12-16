#!/usr/bin/env perl

=head1 NAME

uc_stats.pl

=head1 SYNOPSIS

uc_stats.pl -options

=head1 DESCRIPTION

uc_stats.pl uses UCLUST output (.uc files) as input and does some useful
rearrangements/calculations with it. Three main options (-comp_thresh, -comp_perc, -comp_uniq,
each described in detail under "OPTIONS") are used to determine what uc_stats.pl should do with
the input. A fourth main option (-addtofasta) is a neccessary preliminary step for -comp_perc
and -comp_uniq. Additional tweaks such as outputfile name, table delimiter, print/skip headers,
can be set via sub-options to each main option.

=head1 OPTIONS

=head2 -comp_thresh

Compares clusterings of the same input, but assembled under different UCLUST conditions
(such as different seqIDthresholds or modes).

possible sub-options:

=head3 -i	input filename

=head3 -f	names of additional files that you want to match (see description below)

=head3 -o	set name for output file. if unset, output will be directed to <STDOUT>

=head3 -p	show whether matched read is a centroid (S) or a hit (H) in this clustering (see description below)

=head3 -s	skip headers

=head2

If used on only one .uc file (-f not used), prints out a summarizing table of this file,
showing centroid IDs, cluster no. and cluster sizes.

If used on multiple files (at least one argument in -f), compares clusters of one original .uc file (-i)
to those of additional .uc files (-f). Under these conditions, uc_stats.pl searches each centroid ID
from your main input file (declared under -i) and tries to find it in your additional files. If existing,
prints out the size of the cluster in which this ID was found in, and optionally, whether this ID is a
centroid (S) or a hit (H) in this cluster. 

=head2 -addtofasta

Applies a unique label to each fasta header in a .fasta file.

possible sub-options:

=head3	-i	input filename

=head3	-l	label that should be added to the fasta header (can be any alphanumeric word or character)

=head3	-o	set name for output file. if unset, output will be directed to <STDOUT>

=head2

This is a neccessary preliminary step for -comp_perc and -comp_uniq. These two options will not work
if the clustering was done with unlabeled .fasta files.
Your fasta IDs should start with at least one alphanumeric character BEFORE labeling. 

=head2 -comp_perc

Given one clustering containing reads of two different individuals (pre-labeled with -addtofasta), 
displays the composition of each cluster in percent.

possible sub-options:

=head3	-i	input filename

=head3	-l	labels (must be the same that were used to pre-label .fasta files with -addtofasta)

=head3	-o	set name for output file. if unset, output will be directed to <STDOUT>

=head3	-d	table delimiter (default: "\t")

=head2

Example: You have two .fasta files (A.fasta and B.fasta), each containing reads from one individual.

=head3	1. Pre-label your .fasta files with -addtofasta:
	
		uc_stats.pl -addtofasta -i A.fasta -l A -o A_labeled.fasta
		uc_stats.pl -addtofasta -i B.fasta -l B -o B_labeled.fasta

=head3	2. Concatenate both files:
	
		cat A.fasta B.fasta > all.fasta

=head3	3. Cluster these reads with UCLUST:
	
		usearch -cluster_fast all.fasta <-options>

=head3	4. Check the composition of your clusters with -comp_perc:
	
		uc_stats.pl -comp_perc -i uclust_output.uc -l A B -o uclust_output.comp_perc.csv (-d ",")

=head2 -comp_uniq

Given one clustering containing reads of two or more different individuals (pre-labeled with -addtofasta), 
displays how many clusters are uniquely composed of reads from each single individuals and how many 
clusters share reads from two or more different individuals.

possible sub-options:

=head3	-i	input filename

=head3	-l	labels (must be the same that were used to pre-label .fasta files with -addtofasta)

=head3	-o	set name for output file. if unset, output will be directed to <STDOUT>

=head3	-d	table delimiter (default: "\t")

=head2

Example: You have two .fasta files (A.fasta and B.fasta), each containing reads from one individual.

=head3	1. Pre-label your .fasta files with -addtofasta:
	
		uc_stats.pl -addtofasta -i A.fasta -l A -o A_labeled.fasta
		uc_stats.pl -addtofasta -i B.fasta -l B -o B_labeled.fasta

=head3	2. Concatenate both files:
	
		cat A.fasta B.fasta > all.fasta

=head3	3. Cluster these reads with UCLUST:
	
		usearch -cluster_fast all.fasta <-options>

=head3	4. Check the composition of your clusters with -comp_uniq:
	
		uc_stats.pl -comp_uniq -i uclust_output.uc -l A B -o uclust_output.comp_perc.csv (-d ",")

=cut


use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use List::Util qw(sum);
use Data::Dumper;

my $haaalp;
my $comp_thresh	= 0;
my $addtofasta	= 0;
my @labels;
my $comp_perc	= 0;
my $comp_uniq	= 0;
my $inputfile;
my @files2compare;
my $outfilename;
my $d		= "\t";
my $showsh;
my $skip_headers;


GetOptions (
	"help"			=> \$haaalp,
	"comp_thresh"		=> \$comp_thresh,
	"addtofasta"		=> \$addtofasta,
	"labels=s{,}"		=> \@labels,
	"comp_perc"		=> \$comp_perc,
	"comp_uniq"		=> \$comp_uniq,
	"inputfile=s"		=> \$inputfile,
	"files2compare=s{,}"	=> \@files2compare,
	"outfilename=s"		=> \$outfilename,
	"delimiter=s"		=> \$d,
	"print_SH"		=> \$showsh,
	"skip_header"		=> \$skip_headers,
);

print <<HALP

Usage: uc_stats.pl -<options>

	-h		print this help message

	-comp_thresh	compare  UCLUST output (.uc files) of the same input, but generated with different UCLUST settings.
		-i	input filename
		-f 	additional files to match
		-o	name your output file
		-d	outfile delimiter (default: "\\t")
		-p	show whether matched read is a centroid (S) or a hit (H) in this clustering
			(might be cooler for calculations/plotting)
		-s	skip headers 

	-addtofasta	add a specific label to fasta headers (required for -comp_perc and -comp_uniq)
		-i	input filename
		-l	label that should be added to the fasta header
		-o	name your output file

	-comp_perc	if clustering was done with reads from two different input files (-addtofasta required), 
			gives the composition of each cluster in percent as output.
		-i	input filename
		-l	labels (the same you used to mark your fasta headers in -addtofasta)
		-o	name your output file
		-d	outfile delimiter (default: "\\t")

	-comp_uniq	if clustering was done with reads from two or more different input files (-addtofasta
			required), shows whether clusters are unique for one label or shared as output.
		-i	input filename
		-l	labels (the same you used to mark your fasta headers in -addtofasta)
		-o	name your output file
		-d	outfile delimiter (default: "\\t")

for more detailed instructions, type 'pod2usage -verbose 3 /path/to/uc_stats.pl' into your command line,
or export with any other POD interpreter.

HALP
if $haaalp;
exit if $haaalp;

die <<USAGE


Usage: uc_stats.pl -<options>

	-h		print this help message

	-comp_thresh	compare  UCLUST output (.uc files) of the same input, but generated with different UCLUST settings.
		-i	input filename
		-f 	additional files to match
		-o	name your output file
		-d	outfile delimiter (default: "\\t")
		-p	show whether matched read is a centroid (S) or a hit (H) in this clustering
			(might be cooler for calculations/plotting)
		-s	skip headers 

	-addtofasta	add a specific label to fasta headers (required for -comp_perc and -comp_uniq)
		-i	input filename
		-l	label that should be added to the fasta header
		-o	name your output file

	-comp_perc	if clustering was done with reads from two different input files (-addtofasta required), 
			gives the composition of each cluster in percent as output.
		-i	input filename
		-l	labels (the same you used to mark your fasta headers in -addtofasta)
		-o	name your output file
		-d	outfile delimiter (default: "\\t")

	-comp_uniq	if clustering was done with reads from two or more different input files (-addtofasta
			required), shows whether clusters are unique for one label or shared as output.
		-i	input filename
		-l	labels (the same you used to mark your fasta headers in -addtofasta)
		-o	name your output file
		-d	outfile delimiter (default: "\\t")

for more detailed instructions, type 'pod2usage -verbose 3 /path/to/uc_stats.pl' into your command line,
or export with any other POD interpreter.

USAGE
unless $inputfile;

die "please select one (and only one) of these options: -comp_thresh -addtofasta -comp_perc -comp_uniq\n"
unless ($comp_thresh + $addtofasta + $comp_perc + $comp_uniq) == 1;

#### option -comp_thresh  #####################################################

if ($comp_thresh){

# load centroids, clusternr and clustersize into arrays

	my @centroids;
	my @clusternr;
	my @clustersize;
	my $indexcount = 0;
	
	open my $fh, "<$inputfile" or die $!;
	
	while (<$fh>) {
		my @line = split (/\t/, $_);
		if ($line[0] =~ /C/) {
			$centroids[$indexcount] = $line[8];
			$clusternr[$indexcount] = $line[1];
			$clustersize[$indexcount] = $line[2];
			$indexcount++;
		}
		else {next};
	}
	
	close $fh;

#------------------------------------------------------------------------------
	if (defined($outfilename)){
		open STDOUT, ">$outfilename" or die $!;
	}

	#print $out Dumper(\@centroids);
	#print $out Dumper(\@clusternr);
	#print $out Dumper(\@clustersize);
	
	my %comparehash;
	my %clusterhash;
	my %centroidhash;
	my @matchedline;
	my @matchedcluster;
	my $i;
	my $j;
	
	if (@files2compare){

#### print all the headers-----------------------------------------------------

		unless ($skip_headers){
			print "$inputfile$d$d";
			for my $fileheader (@files2compare) {
				print "$fileheader$d";
				print "$d" if $showsh;
			}
			print "\n";
			print "centroid ID", "$d", "cluster no.", "$d", "size", "$d";

			for my $fileheader2 (@files2compare) {
				print "S/H$d" if $showsh;
				print "size$d";
			}
			print "\n";
		}
		
		for my $file2compare (@files2compare) {
			($clusterhash{$file2compare}, $centroidhash{$file2compare}, $comparehash{$file2compare}) = &loaduc2hash ($file2compare);
		}
#		print $out Dumper(\%clusterhash);
		for my $centroid (@centroids) {
	
#### print stats from inputfile------------------------------------------------
	
			my $printclusternr      = shift @clusternr;
			my $printclustersize    = shift @clustersize;
		        print "$centroid", $d, "$printclusternr", $d, "$printclustersize", "$d";
	
			for my $file2compare (@files2compare){
	
#					print $out Dumper(\%clusterhash);
#					print $out Dumper(\%comparehash);
#					print $out Dumper(\%centroidhash);
				$i = 0;
	
#### match centroids to hits in files2compare----------------------------------
#### and load matched lines into arrays----------------------------------------
	
				foreach (keys $centroidhash{$file2compare}){
					if ($centroidhash{$file2compare}{$_}{line} =~ /$centroid/) {
						$i = 1;
						@matchedline = split(/\t/, $centroidhash{$file2compare}{$_}{line});
						foreach (keys $clusterhash{$file2compare}){
							if ($clusterhash{$file2compare}{$_}{line} =~ /$centroid/) {
								@matchedcluster = split(/\t/, $clusterhash{$file2compare}{$_}{line});
							}
						}
					}
				}
				if ($i == 0){
					$j = 0;
					foreach (keys $comparehash{$file2compare}){
						if ($comparehash{$file2compare}{$_}{line} =~ /$centroid/) {
							$j = 1;
							@matchedline = split(/\t/, $comparehash{$file2compare}{$_}{line});
							foreach (keys $clusterhash{$file2compare}) {
								if ($clusterhash{$file2compare}{$_}{line} =~ /$centroid/) {
									@matchedcluster = split(/\t/, $clusterhash{$file2compare}{$_}{line});
								}
							}
						}
					}
					if ($j == 0){
						@matchedline = "NA";
						@matchedcluster = ("NA", "NA", "NA");
					}
				}
	
#### print stuff out-----------------------------------------------------------
#				print $out Dumper(\@matchedline);
#				print $out Dumper(\@matchedcluster);
				print "$matchedline[0]$d" if $showsh;
				print "$matchedcluster[2]$d";
			}
			print "\n";
		}
	}

#### a little simpler with only one file---------------------------------------

	else {
		print "$inputfile$d$d\n" unless @files2compare;

		while ($indexcount > 0) {
			my $printcentroid 	= shift @centroids;
			my $printclusternr 	= shift @clusternr;
			my $printclustersize 	= shift @clustersize;
			print "$printcentroid", $d, "$printclusternr", $d, "$printclustersize\n";
			$indexcount -= 1;
		}
	}
	close STDOUT if $outfilename;
}

#### option -addtofasta  ######################################################

elsif ($addtofasta){
	die "no label declared!\n" unless @labels;

	open my $fh, "<$inputfile" or die $!;
	if (defined($outfilename)){
		open STDOUT, ">$outfilename" or die $!;
	}
	
	while(<$fh>){
		if ($_ =~ /^>/){
			chomp;
			my @line = split(/\t/, $_);
			$line[0] =~ s/^>/>@labels/;
			print join("\t", @line), "\n";
		}
		else {print};
	}
	close $fh;
	close STDOUT if $outfilename;
}

#### option -comp_perc  #######################################################

elsif ($comp_perc){
	die "no labels declared!\n" unless @labels;

	open my $fh, "<$inputfile" or die $!;

	if (defined($outfilename)){
		open STDOUT, ">$outfilename" or die $!;
	}
	
	my %labels;
	my %input;
	my $sum_labels;
	my $line 		= 0;
	my $cluster	 	= 0;
	
	while (<$fh>){

# dump centroids (S) and hits (H) in a multhash--------------------------------

		if (/^(S|H)(\s+)/){
			$line++;
			if (/\s(\d+)\s/){
				$cluster = $1;
			}
			chomp($input{$cluster}{$line} = $_);
		}
		last if /^C/;
	}
	close $fh;
	
	printf "%s$d%s$d%s\n", "cluster", "sex", "percentage";

# check composition of each cluster--------------------------------------------

	foreach my $key1 (sort {$a <=> $b} keys %input){
		$sum_labels = 0;
		for my $label (@labels){
			$labels{$label} = 0;	
			foreach my $line (values $input{$key1}){
				if ($line =~ /(0|\*)(\s+)(=|\w+|=(\w+)|\*)(\s+)${label}\D/){
					$labels{$label}++;
				}
			}
			$sum_labels += $labels{$label};
		}

# calculate percentage---------------------------------------------------------

		for my $label (@labels){
			printf "%d${d}%s$d%.4f\n", $key1, $label, $labels{$label} / $sum_labels * 100;
		}	
	}
	close STDOUT if $outfilename;
}

#### option -comp_uniq  #######################################################

elsif ($comp_uniq){

	open my $fh, "<$inputfile" or die $!;
	if (defined($outfilename)){
		open STDOUT, ">$outfilename" or die $!;
	}

	my %labels;
	my %input;
	my $line 	= 0;
	my $clusternr 	= 0;
	
	while (<$fh>){
		if (/^(S|H)(\s+)/){
			$line++;
			if (/\s(\d+)\s/){
				$clusternr = $1;
			}
			chomp($input{$clusternr}{$line} = $_);
		}
		last if /^C/;
	}
	close $fh;
	
	#print Dumper(\%input);
	#exit
	
	
	my $sharedcluster;
	my $shared = 0;
	
	printf "%s$d%s\n", "composition", "NumClusters";
	
	for my $label (@labels){
		$labels{$label}{clusterhit} = 0;
	}
	
	foreach my $cluster (sort {$a <=> $b} keys %input){
		$sharedcluster = 0;
		for my $label (@labels){
			$labels{$label}{readhit} = 0;
			foreach my $read (values $input{$cluster}){
				if ($read =~ /(0|\*)(\s+)(=|\w+|=(\w+)|\*)(\s+)${label}\D/){
					$labels{$label}{readhit}++;
					$sharedcluster++;
				}
			last if $labels{$label}{readhit} == 1;
			}
		}
		if ($sharedcluster > 1) {
			$shared++;
		}
		else {
			for my $label (@labels){
				if ($labels{$label}{readhit} == 1){
					$labels{$label}{clusterhit}++;
				}
			}
		}
	}
	
	for my $label (@labels){
		printf "%s$d%d\n", $label, $labels{$label}{clusterhit};
	}
	printf "%s$d%d\n", "shared", $shared;
	close STDOUT if $outfilename;
}

#######################################################################

sub loaduc2hash {
	my $filename = shift @_;
	my %comparefile;
	my %cluster;
	my %centroid;
	my $linecounter 	= 0;
	my $clustercounter	= 0;
	my $centroidcounter 	= 0;

	open my $fh, "<$filename" or die $!;
	while (<$fh>){
		if ($_ =~ /^C/) {
			$cluster{$clustercounter}{line} = $_;
			$clustercounter++;
		}	
		elsif ($_ =~ /^S/) {
			$centroid{$centroidcounter}{line} = $_;
			$centroidcounter++;
		}
		elsif ($_ =~ /^H/) {	
			$comparefile{$linecounter}{line} = $_;
			$linecounter++;
		}
	}
	close $fh;
	return (\%cluster, \%centroid, \%comparefile);
}



