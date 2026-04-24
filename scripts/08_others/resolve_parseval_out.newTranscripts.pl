#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

open IN, $ARGV[1-1] or die $!;

my %new;
my $id;

while (<IN>){
	chomp;
	if (/Unmatched (\S+) transcripts/){
		$id = $1;
		$new{$id} = [];
	}elsif (/\s+\|\s+(\S+)$/ and $id){
		push @{$new{$id}}, $1;
	}elsif (/Locus/ and $id){
		for my $key (sort {$b cmp $a} keys %new){
			print "$key\t".join(";", @{$new{$key}}),"\t",scalar(@{$new{$key}}),"\n";
		}
		$new{$id}= [];
		$id = "";
	}
}
close IN;
