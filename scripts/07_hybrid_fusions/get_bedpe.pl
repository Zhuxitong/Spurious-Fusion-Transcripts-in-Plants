#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

if (@ARGV != 2) {
    die "Usage: perl script.pl <list> <bedpe>\n";
}

my $cistrans = $ARGV[0];
my $stdin = $ARGV[1];

my %list;
open(my $fh_a, '<', $cistrans) or die "Cannot open file '$cistrans': $!\n";
while (<$fh_a>){
	chomp;
	my @l = split/\t/;
	$list{$l[1-1]}++;
}
close $fh_a;


open my $fh_b, $stdin or die $!;
while (<$fh_b>){
	chomp;
	my @l = split/\t/;
	if (exists $list{$l[7-1]}){
		print "$_\n";
	}else{
		next;
	}
}
close $fh_b;
