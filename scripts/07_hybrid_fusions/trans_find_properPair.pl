#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

while (<>){
	chomp;
	my @line = split/\t/;
	my @mh   = @line[3-1..10-1];
	my @zs   = @line[11-1..18-1];
	my $mh_gene1 = $line[3-1];
	my $mh_gene2 = $line[7-1];
	my $mh_gene1_num = $line[4-1];
	my $mh_gene2_num = $line[8-1];
	my $zs_gene1 = $line[11-1];
	my $zs_gene2 = $line[15-1];
	my $zs_gene1_num = $line[12-1];
	my $zs_gene2_num = $line[16-1];

	my $which;
	if ($mh_gene1_num == 1 and $mh_gene2_num == 1){
		$which = 'mh_1_1';
		print $which, "\t", join("\t", @line[1-1,2-1], @mh), "\n";
		#print "--------$_\n";
		next;
	}elsif ($zs_gene1_num == 1 and $zs_gene2_num == 1){
		$which = 'zs_1_1';
		print $which, "\t", join("\t", @line[1-1,2-1], @zs), "\n";
		#print "--------$_\n";
		next;
	}else{
		$which = 'multie';
		print $which,"\t", $_, "\n";
	}
}
