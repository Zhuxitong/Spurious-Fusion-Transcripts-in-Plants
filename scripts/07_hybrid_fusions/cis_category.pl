#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


while (<>){
	chomp;
	my @l = split/\t/;

	# find cis
	if ($l[2-1] =~ 'cis'){
		# find intergenic
		if ($l[3-1] eq "Inter" or $l[7-1] eq "Inter"){
			print "Inter\t$_\n";
		}elsif ($l[3-1] eq $l[7-1] and $l[5-1] eq $l[9-1]){
			print "Same_exon\t$_\n";
		}else {
			print "Diff_exon\t$_\n";
		}
	}else{
		next;
	}
}



