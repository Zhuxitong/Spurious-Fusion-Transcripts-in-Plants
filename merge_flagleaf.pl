#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use List::Util qw/uniq/;

my (%count, %info, %type, %split, %dis, %cov);
while (<>){
	chomp;
	my @line = split/\t/;
	my $id   = join("\t", @line[1-1..4-1]);

	$count{$id}++;
	push @{ $info{$id} }, join("\t", @line[5-1..12-1]);
	push @{ $type{$id} }, $line[13-1];
	$split{ $id } += $line[14-1];
	$dis{ $id } += $line[15-1];
	$cov{ $id } += $line[16-1];
}
# print Dumper \%count;
# print Dumper \%type;

for my $id (keys %count){
	my @info = @{ $info{$id} };
	my @type = @{ $type{$id} };

	my %order_hash = ('high' => 1,'medium' => 2,'low' => 3);
	my @type_sort = sort {$order_hash{$a} cmp $order_hash{$b}} @type;
	print join("\t", $id, $info[1-1], $type_sort[1-1], $split{$id}, $dis{$id}, $cov{$id}),"\n";
	
}
