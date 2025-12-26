#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use List::Util qw/uniq/;


my (%count, %info, %type, %tissue);
while (<>){
	chomp;
	my @line = split/\t/;
	my $id   = join "\t", @line[1-1..4-1];
	$count{$id}++;
	push @{ $info{$id} }, join("\t", @line[5-1..12-1]);
	push @{ $type{$id} }, $line[13-1];
	$tissue{ $id }{$line[17-1]}++;
}

#print Dumper \%count;
# print Dumper \%tissue;

for my $id (keys %count){
	my @info = @{ $info{$id} };
	my %t    = %{ $tissue{$id} };

	my @type = @{ $type{$id} };
	my %order_hash = ('high' => 1,'medium' => 2,'low' => 3);
	my @type_sort = sort {$order_hash{$a} cmp $order_hash{$b}} @type;

	my ($flagleaf, $leaf, $panicle, $root) = (0,0,0,0);
	if (exists $t{flagleaf}){
		$flagleaf = 1;
	}
	if (exists $t{leaf}){
		$leaf = 1;
	}
	if (exists $t{panicle}){
		$panicle = 1;
	}
	if (exists $t{root}){
		$root = 1;
	}
	print join("\t", $id, $info[1-1], $type_sort[1-1], $count{$id}, $flagleaf, $leaf, $panicle, $root), "\n";
	#for my $tissue ('flagleaf', 'leaf', 'panicle', 'root'){
	#	if 
	#}
}
