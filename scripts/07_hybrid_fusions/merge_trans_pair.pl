#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my (%trans, %t);
while (<>){
	chomp;
	my @line = split/\t/;
	
	my $id = $line[1-1];
	my $type = $line[2-1];
	my $pair = join("\t", @line[3-1..10-1]);

	push @{ $trans{$id} }, $pair;
	$t{$id} = $type;
}
#print Dumper \%trans;
#print Dumper \%t;

for my $id (sort keys %trans){
	print join("\t", $id, $t{$id}, @{$trans{$id}}),"\n";
}



#for my $t (keys %t){
#	my @m = @{ $t{$t} };
#	print "$t\t".join("\t", @m),"\n";
#}
