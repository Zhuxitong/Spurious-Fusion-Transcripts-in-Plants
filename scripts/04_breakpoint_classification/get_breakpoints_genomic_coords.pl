#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


while (<>){
	chomp;
	my @a = split/\t/;
	my $id = $a[1-1]."#".$a[2-1];
	my $id1 = $id."_break1";
	my $id2 = $id."_break2";

	print join("\t", $a[3-1], $a[4-1]-1, $a[4-1], $id1, ".", $a[7-1]),"\n";
	print join("\t", $a[5-1], $a[6-1]-1, $a[6-1], $id2, ".", $a[8-1]),"\n";
}
