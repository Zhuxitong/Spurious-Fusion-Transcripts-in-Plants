#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my ($genemap, $in) = @ARGV;


my %geneid;
open IN, $genemap or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	$geneid{$line[2-1]} = $line[1-1];
}
close IN;
#print Dumper \%geneid;


my (%ids, %info, %type);
open IN, $in or die $!;
while (<IN>){
	chomp;
	my @line = split/\t/;
	my $id = (split/_/, $line[4-1])[1-1];
	my $pair_id  = $line[4-1];
	my $pair_num = (split/_/, $line[4-1])[2-1];
	my $type = $line[7-1];

	$ids{$id}{$pair_id}++;
	$type{$id} = $type;
	#push @{ $type{$id} }, $type;

	# process reads with no intersect
	my $gene = 'Unmap';
	if ($line[1-1] eq "ChrUn" and $line[8-1] eq "."){
		$line[16-1] = "ID=Unmap;Parent=Unmap;";
	}
	if ($line[1-1] ne "ChrUn" and $line[8-1] eq "."){
		$line[16-1] = "ID=Inter;Parent=Inter;";
		$gene = 'Inter';
	}
	# get reads info
	my ($exon, $trans) = $line[16-1] =~ /ID=(.*?);Parent=(.*?);/;
	if (exists $geneid{$trans}){
		$gene = $geneid{ $trans };
	}
	push @{ $info{$pair_id}{$gene}} , $exon;
	#print "$_\t$exon\t$trans\n";

}
close IN;
# print Dumper \%info;
#print Dumper \%ids;


for my $id (keys %ids){
	#print "$id\n";
	my $id_1 = $id."_1";
	my $id_2 = $id."_2";
	
	my ($out_gene1, $out_gene2, $out_exon1, $out_exon2) = ('Unmap','Unmap','Unmap','Unmap');
	if (exists $ids{$id}{$id_1}){
		my %gene1 = %{ $info{$id_1} };
		$out_gene1 = join(";", sort keys %gene1);
		my @out;
		for my $x (sort keys %gene1){
			push @out, join(",", @{$gene1{$x}});
		}
		$out_exon1 = join(";", @out);
	}
	if (exists $ids{$id}{$id_2}){
		my %gene2 = %{ $info{$id_2} };
		$out_gene2 = join(";", sort keys %gene2);
		my @out;
		for my $x (sort keys %gene2){
			push @out, join(",", @{ $gene2{$x} });
		}
		$out_exon2 = join(";", @out);
	}
	#print "$id\t$type{$id}\t$out_gene1\t$out_exon1\t$out_gene2\t$out_exon2\n";
	

	my ($gene1_count, $gene2_count, $exon1_count, $exon2_count);
	if ($out_gene1 eq "Unmap" or $out_gene1 =~ /Inter/){
		$gene1_count = 0;
		$exon1_count = 0;
	}else{
		$gene1_count = scalar(split/;/, $out_gene1);
		my @m = split/;/, $out_exon1;
		my @n;
		for my $m (@m){
			push @n, scalar (split/,/, $m);
		}
		$exon1_count = join(";", @n);
	}
	if ($out_gene2 eq "Unmap" or $out_gene2 =~ /Inter/){
		$gene2_count = 0;
		$exon2_count = 0;
	}else{
		$gene2_count = scalar(split/;/, $out_gene2);
		my @m = split/;/, $out_exon2;
		my @n;
		for my $m (@m){
			push @n, scalar (split/,/, $m);
		}
		$exon2_count = join(";", @n);
	}
	print "$id\t$type{$id}\t$out_gene1\t$gene1_count\t$out_exon1\t$exon1_count\t$out_gene2\t$gene2_count\t$out_exon2\t$exon2_count\n";
}
