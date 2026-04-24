#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

open IN, $ARGV[1-1] or die $!;

my (%new,%match);
my %perfect;
my ($id, $gene, $cds_coe, $utr_coe, $overall);
my ($ref_id, $pre_id, $type, $rp);
my ($ref_match, $ref_mismatch, $pre_match, $pre_mismatch);
print "ref_id\tpre_id\tref_exon\tpre_exon\tref_exon_m\tref_exon_mis\tpre_exon_m\tpre_exon_mis\texon_perfect\tref_cds\tpre_cds\tref_cds_m\tref_cds_mis\tpre_cds_m\tpre_cds_mis\tcds_perfect\tref_utr\tpre_utr\tref_utr_m\tref_utr_mis\tpre_utr_m\tpre_utr_mis\tutr_perfect\tgene_perfect\tcds_coe\tutr_coe\toverall_coe\n";

while (<IN>){
	chomp;
	if (/Begin comparison/){
		$id = 'begin';
		%perfect = (CDS => 'No', Exon => 'No', UTR => 'No');
		$gene = 'No';
	}elsif (/reference transcripts:/ and $id){
		my $new_line = <IN>;
		($ref_id) = $new_line =~ /\s+\|\s+(\S+$)/;
	}elsif (/prediction transcripts:/ and $id){
		my $new_line = <IN>;
		($pre_id) = $new_line =~ /\s+\|\s+(\S+$)/;
	}elsif (/(CDS|Exon|UTR) structure comparison/ and $id){
		$type = $1;
		$new{$type} = [];
	}elsif (/(\d+)\s(reference)\s+(CDS|exons|UTR)/ and $id and $type){
		# $3 =~ s/exons/Exon/ if $3 =~ /exons/;
		# push @{$new{$type}}, $2.":".$1;
		push @{$new{$type}}, $1;
		$match{$type}{ref_match} = 'NA';
		$match{$type}{ref_mismatch} = 'NA';
	}elsif (/(\d+)\s(prediction)\s+(CDS|exons|UTR)/ and $id and $type){
		#push @{$new{$type}}, $2.":".$1;
		push @{$new{$type}}, $1;
		$match{$type}{pre_mismatch} = 'NA';
		$match{$type}{pre_match} = 'NA';
	}elsif (/(\d+)\s+match prediction/ and $id and $type) {
		$match{$type}{ref_match} = $1;
	}elsif (/(\d+)\s+don\'t match prediction/ and $id and $type){
		$match{$type}{ref_mismatch} = $1;
	}elsif (/(\d+)\s+match reference/ and $id and $type) {
		$match{$type}{pre_match} = $1;
	}elsif (/(\d+)\s+don\'t match reference/ and $id and $type){
		$match{$type}{pre_mismatch} = $1;
	}elsif (/(CDS|Exon|UTR)\s+structures match perfectly!/ and $id and $type){
		$perfect{$1} = 'perfectly';
	}elsif (/Gene structures match perfectly!/ and $id and $type){
		$gene = 'perfectly';
		($cds_coe, $utr_coe, $overall) = (1,1,1);
	}elsif (/Matching coefficient:\s+(\S+)\s+(\S+)\s+(\S+)/ and $id and $type){
		$cds_coe = $1;
		$utr_coe = $2;
		$overall = $3;
	}elsif (/End comparison/ and $id){
		#print Dumper(\%new),"\n";
		#print Dumper(\%match),"\n";
		#print Dumper(\%perfect),"\n";
		#print join("\t", ($gene, $cds_coe, $utr_coe, $overall)),"\n";
		print "$ref_id\t$pre_id\t";
		for my $key ('Exon','CDS','UTR'){
			print join("\t", @{$new{$key}})."\t".join("\t",($match{$key}{ref_match},$match{$key}{ref_mismatch},$match{$key}{pre_match},$match{$key}{pre_mismatch}))."\t".$perfect{$key}."\t";
		}
		print join("\t", ($gene, $cds_coe, $utr_coe, $overall)),"\n";
		#for my $key (sort {$b cmp $a} keys %new){
		#	print "$key\t".join(";", @{$new{$key}}),"\t",scalar(@{$new{$key}}),"\n";
		#}
		$new{$type}= [];
		$id = "";
	}
}
close IN;
