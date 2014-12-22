#!/usr/bin/env perl


=head1 NAME

compsam.pl - takes one .uc and one .sam file of the same reads as input, calculates how many or which genes the reads of each cluster were mapped to

=head1 SYNOPSIS

compsam.pl -uc <myfile.uc> -sam <myfile.sam> -out <myfile.csv> (-uniq_genes|-num_uniq)

=head1 DESCRIPTION

Output can be given in three forms:

=head3 no specification

Displays the fasta descriptor of the gene that each read of one cluster was mapped to (one entry per read)

=head3 -uniq_genes

Displays the fasta descripter of unique genes that reads of each cluster were mapped to

=head3 -num_uniq

Displays the number of unique genes that reads of each cluster were mapped to

=head2

This script runs horribly slow as it is right now, but it gives you the desired output; just give it some time! ;) However, I might want to re-write this some day when I have the time, maybe include BioPerl Modules like Bio::DC::Sam.

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use Data::Dumper;
use List::MoreUtils qw(uniq);

my $ucfilename;
my $samfilename;
my $outfilename;
my $num_uniq;
my $uniq_genes;

GetOptions(
	"uc=s"		=> \$ucfilename,
	"sam=s"		=> \$samfilename,
	"out=s"		=> \$outfilename,
	"num_uniq"	=> \$num_uniq,
	"uniq_genes"	=> \$uniq_genes,
);


my %ucfilehash 	= loaduc2hash ($ucfilename);
my %samfilehash = loadsam2hash ($samfilename);
my %resulthash;

#print Dumper (\%ucfilehash);
#print Dumper (\%samfilehash);
#exit;

foreach my $cluster (sort {$a <=> $b} keys(%ucfilehash)){
	foreach my $read (values $ucfilehash{$cluster}){
		foreach my $mappedread (values %samfilehash){
			chomp $mappedread;
			if ($mappedread =~ /^$read/){
				my @matchedline = split(/\t/, $mappedread);
				unless (exists($resulthash{$cluster}{$mappedread})){
					$resulthash{$cluster}{$mappedread} = $matchedline[2];
				}
				last;
			}
		}
	}
}

#print Dumper (\%resulthash);
#exit;


my @genes;

if (defined $outfilename){
	open STDOUT, ">$outfilename" or die $!;
}

foreach my $cluster (sort {$a <=> $b} keys(%resulthash)){
	print "c", $cluster;
	foreach my $mappedread (keys $resulthash{$cluster}){
		foreach my $gene ($resulthash{$cluster}{$mappedread}){
			if ($num_uniq || $uniq_genes){
				push(@genes, $resulthash{$cluster}{$mappedread});
			}
			else {
				printf "\t%s", $gene;
			}
		}
	}
	if ($uniq_genes){
		print "\t", join("\t", uniq(@genes));
		print "\n";
	}
	elsif ($num_uniq){
		my $numuniq = uniq(@genes);
		printf "\t%d", $numuniq;
		print "\n";
	}
	else {print "\n"};
}

close STDOUT if $outfilename;


	
###############################################################################

sub loaduc2hash {
	my $filename = shift @_;
	my %filehash;
	my @line;
	my $clusternr;
	my $ucline = 0;

	open my $fh, '<', $filename or die $!;

	while (<$fh>){
		$ucline++;
		my @line = split(/\t/, $_);
		$clusternr = $line[1];
		if ($line[0] =~ /(S|H)/){
			$filehash{$clusternr}{$ucline} = $line[8];
		}
		last if /^C/;
	}
	close $fh;
	return %filehash;
}

###############################################################################

sub loadsam2hash {
	my $filename = shift @_;
	my %filehash;
	my $samline = 0;

	open my $fh, '<', $filename or die $!;

	while (<$fh>){
		next if /^@/;

		$samline++;
		$filehash{$samline} = $_;
	}
	close $fh;
	return %filehash;
}
			

