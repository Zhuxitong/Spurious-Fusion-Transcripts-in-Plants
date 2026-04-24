#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my %gene;
while (<>){
	chomp;
	my @line = split/\t/;
	next if $line[1-1] eq $line[3-1];
	my $id   = $line[1-1]."\t".$line[3-1];
	$gene{$id}++;

	my $id_r = $line[3-1]."\t".$line[1-1];
	if (exists $gene{$id_r}){
		print "---------Found\t$_\n";
	}
}
