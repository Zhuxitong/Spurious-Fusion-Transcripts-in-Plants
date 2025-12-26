#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

while (<STDIN>){
	chomp;
	my @line = split/\t/;
	my @coord1 = split/:/, $line[3-1];
	my @coord2 = split/:/, $line[4-1];
	my $strand1 = (split/\//, $line[5-1])[1-1];
	my $strand2 = (split/\//, $line[6-1])[1-1];
	print join("\t", $coord1[1-1], $coord1[2-1]-1, $coord1[2-1], $coord2[1-1],$coord2[2-1]-1,$coord2[2-1],$line[1-1]."_".$line[2-1], ".", $strand1, $strand2, @line[7-1..18-1]),"\n";
}
