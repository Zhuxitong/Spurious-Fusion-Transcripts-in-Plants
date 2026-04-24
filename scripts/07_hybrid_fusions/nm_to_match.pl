#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

while (<>){
	chomp;
	my @line = split/\t/;
	my @out;
	for my $nm (@line[2-1..5-1]){
		if ($nm eq 'NA'){
			$nm = 1;
		}
		if ($nm == 0){
			push @out, 'Y';
		}else{
			push @out, 'N';
		}
	}
	print join("\t", $line[1-1], @out),"\n";
}
