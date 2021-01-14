#!/usr/bin/perl -w

@gene = split ",",$ARGV[0];
$result = $ARGV[1]; ## output directory
$ref_fa = $ARGV[2];
$ref_annot = $ARGV[3];
$broad_region = 0;
##
open CO,">$result/combine.fa" or die $!;
foreach $gene (@gene){
	## get gene bed
	$txt =  `grep -w $gene $ref_annot | cut -f 1,4,5,6,7 |sed 's/chr//g' | head -n 1`;chomp($txt);
	# chr7    72300004        72307280        .       +
	open O1,">$result/$gene.bed";
	print O1 "#chr\tstart\tend\tID\tINFO\tstrand\n";
	@tmp = split "\n",$txt;
	foreach $txt (@tmp){
		@txt = split "\t",$txt;
		$s = $txt[1]-$broad_region; 
		$e = $txt[2]+$broad_region; 
		$id= $gene;
		$chr = $txt[0];
		$strand = $txt[4];
		if($s<0){$s=1;}
		print O1 "$chr\t$s\t$e\t$id\t.\t$strand\n";
	}
	close O1;
	$cmd = "seqtk subseq $ref_fa $result/$gene.bed";
	$txt = `$cmd | tail -n 1`; chomp($txt);
	$strand = `cut -f 6 $result/$gene.bed | tail -n 1`;chomp($strand);
	if($strand eq "-"){
		$txt =  uc(reverse($txt));
		$txt =~ tr/ATCG/TAGC/;
	}
	print CO ">$gene\n$txt\n";
}
close CO;

