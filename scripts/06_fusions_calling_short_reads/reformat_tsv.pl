#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

while (<>){
	chomp;
	my @l = split/\t/;
	next if /^#gene/;

	my $gene1 = $l[1-1];
	my $gene2 = $l[2-1];
	my $break1 = $l[5-1];
	my $break2 = $l[6-1];
	my $strand1 = $l[3-1];
	my $strand2 = $l[4-1];

	my $site = 'inter';
	my $chr1 = (split/:/, $break1)[1-1];
	my $chr2 = (split/:/, $break2)[1-1];
	if ($chr1 eq $chr2){
		$site = 'intra';
	}

	my $loc1  = (split/\//, $l[7-1])[1-1];
	my $loc2  = (split/\//, $l[8-1])[1-1];

	my ($splice1, $splice2) = ('N', 'N');
	$splice1 = 'S' if $l[7-1] =~ /splice-site/;
	$splice2 = 'S' if $l[8-1] =~ /splice-site/;
	my $splice = $splice1.'/'.$splice2;

	my $frame = $l[16-1];

	my $adj = '';
	my $shs = '';

	my $type = $l[9-1];
	my $conf = $l[15-1];

	my $split = $l[10-1]+$l[11-1];
	my $dis   = $l[12-1];
	my $cov   = $l[13-1] + $l[14-1];

	# print join("\t", $gene1, $gene2, $break1, $break2, $strand1, $strand2, $site, $splice, $loc1, $loc2, $l[7-1], $l[8-1], $type, $conf, $split, $dis, $cov),"\n";
	print join("\t", $gene1, $gene2, $break1, $break2, $strand1, $strand2, $site, $splice, $loc1, $loc2, $type, $conf, $frame, $split, $dis, $cov),"\n";
}
