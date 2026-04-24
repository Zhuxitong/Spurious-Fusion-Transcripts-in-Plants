#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

while (<>){
	chomp;
	my @line = split/\t/;
	my $mh_gene1 = $line[3-1];
	my $mh_gene2 = $line[7-1];
	my $mh_gene1_num = $line[4-1];
	my $mh_gene2_num = $line[8-1];
	my $zs_gene1 = $line[11-1];
	my $zs_gene2 = $line[15-1];
	my $zs_gene1_num = $line[12-1];
	my $zs_gene2_num = $line[16-1];


	my $type = 'Pass';
	# unmap or intergenic
	if ($line[2-1] eq "MHZStrans"){
		if ($mh_gene1 =~ /Unmap|Inter/ or $zs_gene2 =~ /Unmap|Inter/){
			$type = 'Intergenic';
		}elsif ($mh_gene1_num > 1 and $zs_gene2_num > 1){
			$type = 'Pass_multipleGenes';
		}
	}elsif ($line[2-1] eq "ZSMHtrans"){
		if ($mh_gene2 =~ /Unmap|Inter/ or $zs_gene1 =~ /Unmap|Inter/){
			$type = 'Intergenic';
		}elsif ($mh_gene2_num >1 and $zs_gene1_num >1){
			$type = 'Pass_multipleGenes';
		}
	}else{
		print "Unknown type\t$_\n";
	}
	print "$type\t$_\n";
}
