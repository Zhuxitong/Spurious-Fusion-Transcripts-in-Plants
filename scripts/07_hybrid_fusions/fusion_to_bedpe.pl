#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


my ($bed_fh, $fusion_fh) = @ARGV;
my (%bed, %s);
open IN, $bed_fh or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	$bed{ $line[4-1] } = join("\t", @line[1-1..3-1]);
	$s{ $line[4-1] } = $line[6-1];
}
close IN;
# print Dumper \%s;

open IN, $fusion_fh or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	my $g1 = $line[1-1];
	my $g2 = $line[2-1];
	if (exists $bed{$g1} and exists $bed{$g2}){
		# for ture fusions
		print join("\t", $bed{$g1}, $bed{$g2}, $line[3-1], '.', $s{$g1}, $s{$g2}),"\n";
		# for random. the random id
	}else{
		print "Error\n";
	}
}
close IN;
