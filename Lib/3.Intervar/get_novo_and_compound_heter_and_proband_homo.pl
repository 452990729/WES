#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script designed for retriving denovo, compound heterozygous and parents heter/proband homo mutation (MAF < 0.01).
        Author: zhoujj2013\@gmail.com 
        Usage: $0 <proband_col> <father_col> <mother_col> <sample_count> <AllIntergrateMutation.txt> <prefix>
        Example: perl get_novo_and_compound_heter_and_proband_homo.pl 1 3 2 3 AllIntergrateMutation.txt test

USAGE
print "$usage";
exit(1);
};

my $proband_col = shift;
my $father_col = shift;
my $mother_col = shift;

$proband_col = $proband_col -1;
$father_col = $father_col - 1;
$mother_col = $mother_col - 1;

my $sample_count = shift;

my $all_vcf_f = shift;

my $prefix = shift;

my %proband;
my %father;
my %mother;

open OUT1,">","./$prefix.denovo.mut.txt" || die $!;
open OUT2,">","./$prefix.parent.heter.proband.homo.mut.txt" || die $!;
open IN,"$all_vcf_f" || die $!;
my $header = <IN>;
print OUT1 "$header";
print OUT2 "$header";
while(<IN>){
	chomp;
	my @t = split /\t/;

	my @p = split /:/,$t[$proband_col];
	my @f = split /:/,$t[$father_col];
	my @m = split /:/,$t[$mother_col];
	
	# filter by MAF < 0.01	
	my $exac = $t[31+$sample_count];
	my $esp6500 = $t[40+$sample_count];
	my $g10002015 = $t[41+$sample_count];

	my $flag = 0;
	my @maf = ($exac,$esp6500,$g10002015);
	if($exac eq "." && $esp6500 eq "." && $g10002015 eq "."){
		$flag++;
	}else{
		# if MAF > 0.01 in one dataset, filter the mutations.
		foreach my $maf (@maf){
			next if($maf eq ".");
			if($maf > 0.01){
				$flag = 0;
				last;
			}else{
				$flag++;
			}
		}
	}
	next if($flag == 0);
	
	#print "$exac\t$esp6500\t$g10002015\n";
	
	# find denovo mutations
	my $p_gt = $p[0];
	my $f_gt = $f[0];
	my $m_gt = $m[0];

	my ($p_ref_c,$p_alt_c) = split /,/,$p[1];

	# only select denovo heterozygous mutation
	if(($p_gt ne "0/0" && $p_gt ne "1/1" && $p_gt ne "1|1") && ($f_gt eq "0/0" && $m_gt eq "0/0")){
		if($p[2] >= 10 && $p_alt_c >= 5){
			print OUT1 join "\t",@t;
			print OUT1 "\n";
		}
	}

	# find father heter, mother heter, but proband homo
	if($p_gt eq "1/1" && ($f_gt eq "0/1" && $m_gt eq "0/1")){
		#print "$t[$proband_col]\t$t[$father_col]\t$t[$mother_col]\n";
		print OUT2 join "\t",@t;
		print OUT2 "\n";
	}

	# store proband for compound heterozygous analysis
	my $chr = $t[0+$sample_count];
	my $start = $t[1+$sample_count];
	my $end = $t[2+$sample_count];
	my $ref = $t[3+$sample_count];
	my $alt = $t[4+$sample_count];
	my $geneid = $t[5+$sample_count];

	my $snv_id = "$chr:$start-$end:$ref>$alt";
	if($p_gt eq "0/1"){
		$proband{$geneid}{$snv_id} = 1;
	}

	if($f_gt eq "0/1"){
		$father{$geneid}{$snv_id} = 1;
	}

	if($m_gt eq "0/1"){
		$mother{$geneid}{$snv_id} = 1;
	}
}
close IN;
close OUT1;
close OUT2;


my %compound_heter;

foreach my $gid (keys %proband){
	my @f_snvs;
	my @m_snvs;
	foreach my $snv_id (keys %{$proband{$gid}}){
		if(exists $father{$gid}{$snv_id} && !(exists $mother{$gid}{$snv_id})){
			push @f_snvs,$snv_id;
			next;
		}elsif(exists $mother{$gid}{$snv_id} && !(exists $father{$gid}{$snv_id})){
			push @m_snvs,$snv_id;
			next;
		}
	}
	if(scalar(@f_snvs) > 1 && scalar(@m_snvs)){
		foreach my $mut_id (@f_snvs){
			$mut_id = $mut_id.":$gid";
			$compound_heter{$mut_id} = 1;
		}
		foreach my $mut_id (@m_snvs){
			$mut_id = $mut_id.":$gid";
			$compound_heter{$mut_id} = 1;
		}
	}
}


open OUT3,">","./$prefix.compound_heter.txt" || die $!;
open IN,"$all_vcf_f" || die $!;
$header = <IN>;
print OUT3 "$header";
while(<IN>){
	chomp;
	my @t = split /\t/;
	
	# generate snv_id 
	my $chr = $t[0+$sample_count];
	my $start = $t[1+$sample_count];
	my $end = $t[2+$sample_count];
	my $ref = $t[3+$sample_count];
	my $alt = $t[4+$sample_count];
	my $geneid = $t[5+$sample_count];

	my $snv_id = "$chr:$start-$end:$ref>$alt:$geneid";

	if(exists $compound_heter{$snv_id}){
		print OUT3 join "\t",@t;
		print OUT3 "\n";
	}
}
close IN;
close OUT3;
