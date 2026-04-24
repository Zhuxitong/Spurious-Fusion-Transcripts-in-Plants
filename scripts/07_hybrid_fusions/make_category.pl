#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my %c;
while (<>){
	chomp;
	my @line = split/\t/;
	$c{ join("\t", @line[2-1..5-1])}++
}
#print Dumper \%c;


for my $c (sort keys %c){
	print "$c\t$c{$c}\n";
}
