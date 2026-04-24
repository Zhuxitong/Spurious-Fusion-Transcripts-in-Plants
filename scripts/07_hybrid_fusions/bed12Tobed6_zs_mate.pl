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
	my $fake_coor = "ChrUn\t0\t1\t";
	my $fake_score= "\t.\t+";

	if ($cis{$id} eq "ZScis"){
		print "$r1\t$cis{$id}\n$r2\t$cis{$id}\n";
	}elsif ($cis{$id} eq "ZSMHtrans"){
		my $out;
		if ($l[4-1] eq "."){
			$out = $fake_coor.$id."_2".$fake_score;
		}else{
			$out = $r2;
	}
		print "$r1\t$cis{$id}\n$out\t$cis{$id}\n";
	}elsif ($cis{$id} eq "MHZStrans"){
		my $out;
		if ($l[1-1] eq "."){
			$out = $fake_coor.$id."_1".$fake_score;
		}else{
		    $out = $r1;
		}
		print "$out\t$cis{$id}\n$r2\t$cis{$id}\n";
	}
}
close IN;
