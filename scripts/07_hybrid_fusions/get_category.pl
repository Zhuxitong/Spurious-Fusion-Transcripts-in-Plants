#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

while (<>){
	chomp;
	my @line = split/\t/;
	my $id   = $line[1-1];
	my $t    = join("", @line[2-1..5-1]);
	if ($t eq "YYNN"){
		print "$id\tMHcis\n";
	}elsif ($t eq "NNYY"){
		print "$id\tZScis\n";
	}elsif ($t eq "YNNY"){
		print "$id\tMHZStrans\n";
	}elsif ($t eq "NYYN"){
		print "$id\tZSMHtrans\n";
	}else{
		next;
	}
}

