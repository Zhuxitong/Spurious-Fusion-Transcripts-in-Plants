#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


my ($trans_file, $cis_file) = @ARGV;

my %pair;

open IN, $trans_file or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	my $id   = $line[1-1]."\t".$line[2-1];
	my $trans = $line[4-1];
	my $trans_type = $line[3-1];

	$pair{$id}{trans} = $trans;
	$pair{$id}{trans_type} = $trans_type;
}
close IN;
#print Dumper \%pair;


open IN, $cis_file or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	my $id = $line[1-1]."\t".$line[2-1];
	my $cis = $line[4-1];
	my $cis_type = $line[3-1];

	$pair{$id}{cis} = $cis;
	$pair{$id}{cis_type} = $cis_type;
}
close IN;
# print Dumper \%pair;


for my $id (keys %pair){
	my ($trans, $cis) = (0,0);
	my ($trans_type, $cis_type) = ('NA','NA');

	if (exists $pair{$id}{trans}){
		$trans = $pair{$id}{trans};
	}
	if (exists $pair{$id}{cis}){
		$cis = $pair{$id}{cis};
	}
	if (exists $pair{$id}{trans_type}){
		$trans_type = $pair{$id}{trans_type};
	}
	if (exists $pair{$id}{cis_type}){
		$cis_type = $pair{$id}{cis_type};
	}
	print $id,"\t",join("\t", $trans, $cis, $trans_type, $cis_type),"\n";
}
