#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use List::Util qw/uniq/;

while (<>){
	chomp;
	my @line = split/\t/;
	
	# filter 1 repeat only
	next if $line[18-1] < 2;

	# reorder rank;
	my @rank = uniq(split/;/, $line[13-1]);
	my %order_hash = ('high' => 1,'medium' => 2,'low' => 3);
	my @rank_sort = sort {$order_hash{$a} cmp $order_hash{$b}} @rank;
	# print "$_\n";
	#print Dumper \%order_hash;
	#print "@rank\t||\t@rank_sort\n";
	$line[13-1] = $rank_sort[1-1];

	print join("\t", @line[1-1..16-1]),"\n";
}
