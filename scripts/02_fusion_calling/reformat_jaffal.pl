#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;


my (%s, %r, %span_read, %strand, %add);
while (<>){
	chomp;
	my @l = split/,/;

	my @genes = split/:/, $l[2-1];
	my ($break1, $break2) = ($l[3-1].":".$l[4-1], $l[6-1].":".$l[7-1]);;
	my ($s1, $s2) =($l[5-1], $l[8-1]);
	my $sample = $l[1-1];
	my $span_read = $l[11-1];
	my $reads  = $l[15-1];
	my $type = $l[17-1];

	my $chr  = 'Intra';
	if ($l[9-1] eq "Inf"){
		$chr = "Inter";
	}
	my $info = $l[13-1]."\t".$l[12-1]."\t".$chr;

	my $id = join("\t", @genes, $break1, $break2, $type);
	push @{ $s{$id} }, $sample;
	push @{ $r{$id} }, $reads;
	$span_read{ $id } += $span_read;
	$strand{ $id } = $s1."\t".$s2;
	$add{ $id} = $info;
}

#print Dumper \%s;
#print Dumper \%r;
#print Dumper \%span_read;
#print Dumper \%strand;


for my $id (keys %s){
	my @sample = @{ $s{$id} };
	my @reads  = @{ $r{$id} };
	print join("\t", $id, $strand{$id}, $span_read{$id}, scalar @sample, join(";", @sample), join(";",@reads ), $add{$id}),"\n";
}
