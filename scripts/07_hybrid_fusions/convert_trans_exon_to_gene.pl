#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use List::Util qw(uniq);

my (%pair, %t);
while (<>){
	chomp;
	my @line = split/\t/;
	my $pair = join("\t", sort($line[1-1], $line[3-1]));
	# my $pair = $line[1-1]."\t".$line[3-1];
	my $type = $line[5-1];
	my $num  = $line[6-1];
	
	$pair{$pair} += $num;
	push @{ $t{$pair} }, $type;
}
#print Dumper \%pair;
#print Dumper \%t;


for my $pair (keys %pair){
	my @t = @{ $t{$pair} };
	my @t_uniq = uniq(@t);
	print $pair,"\t",join(";", sort @t_uniq), "\t", $pair{$pair}, "\n";
}
