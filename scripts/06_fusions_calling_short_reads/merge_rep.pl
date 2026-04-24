#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


my (%info, %count, %con,%split, %dis, %cov, %sample);
while (<>){
	chomp;
	my @line = split/\t/;

	my $id = join("\t", @line[1-1..4-1]);
	push @{ $info{ $id } }, join("\t", @line[5-1..11-1,13-1]);
	$count{ $id } ++;
	push @{ $con{ $id } }, $line[12-1];
	$split{$id} += $line[14-1];
	$dis{$id} += $line[15-1];
	$cov{$id} += $line[16-1];
	push @{ $sample{$id} }, $line[17-1];
}

# print Dumper \%count;


for my $id (keys %count){
	my @info = @{ $info{$id} };
	my $con  = join ";", @{ $con{$id} };
	my $sample = join ";", @{ $sample{$id} };
	print "$id\t",join("\t", $info[1-1], $con, $split{$id}, $dis{$id}, $cov{$id}, $sample, $count{$id}),"\n";
	#print "\n";

=pod
	next if scalar @info == 1;
	my $c = 1;
	for my $gene (@info){
		print "$c\t$id\t$gene\n";
		$c++;
	}
	print "\n";
=cut
}
