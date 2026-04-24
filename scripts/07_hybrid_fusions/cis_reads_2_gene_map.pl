#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


while (<>){
	chomp;
	my @line = split/\t/;

	my @gene1 = split/;/, $line[4-1];
	my @gene2 = split/;/, $line[8-1];
	my @exon1 = split/;/, $line[6-1];
	my @exon2 = split/;/, $line[10-1];
	
	my (@mate1, @mate2);
	for (my $x=0; $x<scalar(@gene1); $x++){
		my @exon = split/,/, $exon1[$x];
		for my $exon (@exon){
			#push @mate1, $gene1."\t".$exon;
			print $gene1[$x]."\t".$exon."\n";
		}
	}
	for (my $y=0; $y<scalar(@gene2); $y++) {
		my @exon = split/,/, $exon2[$y];
		for my $exon (@exon){
			#push @mate2, $gene2."\t".$exon;
			print $gene2[$y]."\t".$exon."\n";
		}
	}
#	print Dumper \@mate1;
#	print Dumper \@mate2;
}
