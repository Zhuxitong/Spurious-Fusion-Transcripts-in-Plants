#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my %pair;
while (<>){
	chomp;
	my @line = split/\t/;

	next if ($line[5-1] > 1 or $line[9-1] > 1);

	my $gene1 = $line[4-1];
	my $gene2 = $line[8-1];
	my @exon1 = split/,/, $line[6-1];
	my @exon2 = split/,/, $line[10-1];

	my (@mate1, @mate2);
	for my $exon1 (@exon1){
		push @mate1, $gene1."\t".$exon1;
	}
	for my $exon2 (@exon2){
		push @mate2, $gene2."\t".$exon2;
	}

	foreach my $x (@mate1){
		foreach my $y (@mate2){
			push @{ $pair{ $x."\t".$y."\t".$line[3-1] } }, $line[2-1];
		}
	}

	#print Dumper \@mate1;
	#print Dumper \@mate2;
	#print "------------------\n";	
}
#print Dumper \%pair;

for my $pair (keys %pair){
	my @reads = @{ $pair{$pair} };
	print $pair,"\t",scalar(@reads),"\t", join(";", @reads),"\n";
}
