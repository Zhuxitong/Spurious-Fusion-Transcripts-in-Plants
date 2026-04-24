#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my ($hybird_file, $control_file) = @ARGV;

my %pair;

open IN, $hybird_file or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	my $id   = $line[1-1]."\t".$line[2-1];
	my $h_trans = $line[3-1];
	my $h_cis  = $line[4-1];

	$pair{$id}{h_trans} = $h_trans;
	$pair{$id}{h_cis}   = $h_cis;
}
close IN;


open IN, $control_file or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	my $id = $line[1-1]."\t".$line[2-1];
	my $c_trans = $line[3-1];
	my $c_cis = $line[4-1];

	$pair{$id}{c_trans} = $c_trans;
	$pair{$id}{c_cis}  = $c_cis;
}
close IN;
# print Dumper \%pair;


for my $id (keys %pair){
        my ($h_trans, $h_cis, $c_trans, $c_cis) = (0,0,0,0);

        if (exists $pair{$id}{h_trans}){
                $h_trans = $pair{$id}{h_trans};
        }
        if (exists $pair{$id}{h_cis}){
                $h_cis = $pair{$id}{h_cis};
        }
        if (exists $pair{$id}{c_trans}){
                $c_trans = $pair{$id}{c_trans};
        }
        if (exists $pair{$id}{c_cis}){
                $c_cis = $pair{$id}{c_cis};
        }
		print $id,"\t",join("\t", $h_trans, $c_trans, $h_cis, $c_cis),"\n";
}
