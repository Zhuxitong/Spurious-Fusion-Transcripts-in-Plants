#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


my ($gene_file, $pair_file) = @ARGV;
my %chr;
open IN, $gene_file or die $!;
while (<IN>) {
	chomp;
	my @line = split/\t/;
	push @{ $chr{$line[1-1]} }, $line[4-1];
}
close IN;
# print Dumper \%chr;


open IN, $pair_file or die $!;
while (<IN>){
	chomp;
	my @line  = split/\t/;
	my ($chr1) = $line[1-1] =~ /_(\d{2})G/;
	my ($chr2) = $line[2-1] =~ /_(\d{2})G/;

	my @g_chr1 = @{ $chr{'Chr'.$chr1} };
	my @g_chr2 = @{ $chr{'Chr'.$chr2} };

	my $index1 = int( rand(scalar @g_chr1) );
	my $index2 = int( rand(scalar @g_chr2) );

	# print join("\t", $_, $chr1, $chr2, $g_chr1[$index1], $g_chr2[$index2]),"\n";
	print "$g_chr1[$index1]\t$g_chr2[$index2]\n";
}
close IN;
