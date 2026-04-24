#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

open IN, $ARGV[1-1] or die $!;

my (%ref, %pre);
my ($id, @r, @p);
my $current_key;

while (<IN>){
	chomp;
	if (/Locus/){
		$_ =~ /\|---- Locus: seqid=(\S+) range=(\S+)/;
		print "$1:$2\t";
	}elsif (/(reference genes|prediction genes)/){
		$current_key = $1;
		$ref{$current_key} = [];
	}elsif (/^\|\s+(\S+)/ and $current_key){
		push @{$ref{$current_key}}, $1;
	}elsif (/\|----------$/){
		#print Dumper(\%ref);
		my @out;
		for my $key (sort {$b cmp $a}  keys %ref){
			# push @out, "$key:".join(";", @{$ref{$key}});
			push @out, join(";", @{$ref{$key}});
			my $c = 0;
			for my $n (@{$ref{$key}}){
				$c++ if $n !~ /None/;
			}
			push @out, $c;
		}
		print join("\t",@out),"\n";
	}
}
#print Dumper(\%ref);
close IN;
