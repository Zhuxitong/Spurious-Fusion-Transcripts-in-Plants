#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my (%id, %cov, %nm, %map, %final);

while (<>){
	next if /^@/;

	chomp;
	my @line = split/\t/;

	# cal coverage
	my $cov='NA';
	if ($line[6-1] !~ /\*/){
		#print "$line[6-1]\n";
		my $seq_len = length($line[10-1]);
		my $unmap   = 0;
		my @cigar_op = split/\d+/, $line[6-1];
		shift @cigar_op;
		my @cigar_nu = split/\D+/, $line[6-1];
		for (my $i=0; $i<scalar(@cigar_op); $i++){
			if ($cigar_op[$i] eq 'S' or $cigar_op[$i] eq 'H'){
				$unmap = $unmap + $cigar_nu[$i];
			}else{
				next;
			}
		}
		$cov = ($seq_len - $unmap) / $seq_len;
		#print Dumper \@cigar_nu;
		#print "$cov\n";
	}
	
	# get NM
	my $nm = 'NA';
	for my $x (@line){
		if ($x =~ /NM:i/){
			$nm = (split/:/, $x)[3-1];
		}
	}
	#print "$nm\t$_\n";

	my $flag = $line[2-1];
	
	# judge if read is mapped or not
	my $align = 'Y';
	if ($flag & 0x4){
		$align = 'N';
	}
	
	# get final matrix
	my $final = 'Y';
	if ($cov eq 'NA' or $cov<0.9 or $nm eq 'NA' or $nm>0 or $align eq 'N'){
		$final = 'N';
	}

	# pair info
	if ($flag & 0x40){
		#print "first\t$_\n";
		push @{ $id{$line[1-1]} }, $line[1-1]."_1";
		$cov{ $line[1-1]."_1" } = $cov;
		$nm{ $line[1-1]."_1" } = $nm;
		$map{ $line[1-1]."_1" } = $align;
		$final{ $line[1-1]."_1"} = $final;
	}elsif ($flag & 0x80){
		push @{ $id{$line[1-1]} }, $line[1-1]."_2";
		$cov{ $line[1-1]."_2" } = $cov;
		$nm{ $line[1-1]."_2" } = $nm;
		$map{ $line[1-1]."_2" } = $align;
		$final{ $line[1-1]."_2"} = $final;
	}
}
# print Dumper \%id;


for my $id (sort keys %id){
	my @mate = @{ $id{$id} };
	print "$id\t";
	my @out;
	for my $i (sort @mate){
		#print join("\t", $i, $map{$i}, $cov{$i}, $nm{$i});
		push @out, ($map{$i}, $cov{$i}, $nm{$i}, $final{$i});
	}
	print join("\t", @out), "\n";
}
