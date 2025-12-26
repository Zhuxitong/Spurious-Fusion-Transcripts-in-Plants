#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


my ($t_fh, $jaffal_fh) = @ARGV;

my %pc;
open IN, $t_fh or die $!;
while (<IN>){
	chomp;
	my @l = split/\t/;
	my @sample = split/;/, $l[10-1];
	my @reads  = split/;/, $l[11-1];
	my $id = $sample[1-1]."|".$l[1-1].":".$l[2-1]."|".$reads[1-1];
	$pc{$id} = $l[5-1];
}
close IN;
# print Dumper \%pc;


open IN, $jaffal_fh or die $!;
while (<IN>){

	chomp;
	my @l = split/,/;

	my $genes = $l[2-1];
    # my ($break1, $break2) = ($l[3-1].":".$l[4-1], $l[6-1].":".$l[7-1]);;
    # my ($s1, $s2) =($l[5-1], $l[8-1]);
    my $sample = $l[1-1];
    # my $span_read = $l[11-1];
    my $reads  = $l[15-1];
	# my $type = $l[17-1];
	my $id = $sample."|".$genes."|".$reads;
	my $break = $l[16-1];

	if (exists $pc{$id}){
		print join("\t", $id, $break-1, $break+1, $pc{$id}),"\n";
	}
}
close IN;
