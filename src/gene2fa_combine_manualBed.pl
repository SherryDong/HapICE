#!/usr/bin/perl -w

@gene = split ",",$ARGV[0];
@input_bed = split ",",$ARGV[1];
$result = $ARGV[2]; ## output directory
$ref_fa = $ARGV[3];
##
open CO,">$result/combine.fa" or die $!;
foreach $i (0..$#gene){
	$gene = $gene[$i];
	$bed = $input_bed[$i];
	$out_bed = "$result/$gene.bed";
	if($bed ne $out_bed){
		system "cp $bed $out_bed";
	}
	$cmd = "seqtk subseq $ref_fa $out_bed";
	$txt = `$cmd | tail -n 1`; chomp($txt);
	$strand = `cut -f 6 $out_bed | tail -n 1`;chomp($strand);
	if($strand eq "-"){
		$txt =  uc(reverse($txt));
		$txt =~ tr/ATCG/TAGC/;
	}
	print CO ">$gene\n$txt\n";
}
close CO;

