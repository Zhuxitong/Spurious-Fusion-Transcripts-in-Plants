#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my ($cis_fh ,$bedpe) = @ARGV;

open IN, $cis_fh or die $!;
my %cis;
while (<IN>){
	chomp;
	my @line = split/\t/;
	$cis{ $line[1-1] } = $line[2-1];
}
close IN;
#print Dumper \%cis;


open IN, $bedpe or die $!;
while (<IN>){
	chomp;
	my @l = split/\t/;
	my $id = $l[7-1];
	my $r1 = join("\t", @l[1-1..3-1],$l[7-1]."_1",".",$l[9-1]);
	my $r2 = join("\t", @l[4-1..6-1],$l[7-1]."_2",".",$l[10-1]);
	if ($cis{$id} eq "MHcis"){
		print "$r1\t$cis{$id}\n$r2\t$cis{$id}\n";
	}elsif ($cis{$id} eq "MHZStrans"){
		print "$r1\t$cis{$id}\n";
	}elsif ($cis{$id} eq "ZSMHtrans"){
		print "$r2\t$cis{$id}\n";
	}
}
close IN;

