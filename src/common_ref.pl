#!/usr/bin/perl -w

$aln = $ARGV[0];
$gene = $ARGV[1];
$each = $ARGV[2]+2; ## gene number + match + empty line
$gene_bed = $ARGV[3];
$out1 = $ARGV[4];
$out2 = $ARGV[5];
$out3 = $ARGV[6];
$tmp = `cat $gene_bed | tail -n 1`;chomp($tmp);
@info = split "\t",$tmp;
## 
#perl src/common_ref.pl $combine_ref_dir/combine.aln $true_gene 2 $combine_ref_dir/$true_gene.bed $combine_ref_dir/combine_mod.info $combine_ref_dir/combine_mod.fa $combine_ref_dir/combine_mod_region.bed
open F,$aln or die $!;
$line_num = 0;
$all_seq = "";
$all_mark = "";
while(<F>){
	chomp;
	$line = $_;
	$line_num++;
	if($line_num<3){next;}
	$line_e = ($line_num-3)%$each;
	if($line_e == $each-1){ ## last one
		if($line =~ /^                (.*)$/){
			$mark = $1;
			$all_mark .= $mark;
		}
	}else{
		if($line_e < $each-1){ ## first one ## modify
			if($line =~ /^$gene\s+([A-Za-z-]+)/){
				$seq = $1;
				$all_seq .= $seq;
			}elsif($line =~ /^(.*)\s+([A-Za-z-]+)/){
				$other_gene = $1;
				$seq = $2;
				$other_gene =~ s/\s//g;
				$all_other_seq{$other_gene} .= $seq;
			}	
		}
	}
}
close F;
@all_seq  = split "",$all_seq;
@all_mark = split "",$all_mark;
@all_other_gene = keys %all_other_seq; $n = join "\t",@all_other_gene;
foreach $other_gene(@all_other_gene){
	@{$other_gene} = split "",$all_other_seq{$other_gene};		
}
#print $#all_seq."\t$#all_mark\n";
################
open O1,">$out1";
open O2,">$out2";
open O3,">$out3"; ## 0-based
print O1 "#combine_pos\tori_pos\tgenome_pos\t$gene\t$n\tmark\n";
$ori_i = 0;
$strand = $info[5];
if($strand eq "-"){$start_pos = $info[2];}else{$start_pos=$info[1];}
$mod_start = 0;
$mod_end   = 0;
foreach $i (0..$#all_seq){
	undef(@n);
	foreach $other_gene(@all_other_gene){
		push(@n,@{$other_gene}[$i]);
	}
	$n = join "\t",@n;
	if($all_seq[$i] ne '-'){$ori_i++;}
	if($strand eq "-"){$real_pos = $start_pos-$ori_i+1;}else{$real_pos = $start_pos+$ori_i;}
	print O1 $i."\t$ori_i\t$real_pos\t".$all_seq[$i]."\t$n\t".$all_mark[$i]."\n";
	if($all_mark[$i] eq "*"){
		$final_seq .= $all_seq[$i];
	}else{
		$final_seq .= "N";
	}
	if($ori_i == 1){$mod_start = $i-1;next;}
#print $ori_i."\t".$mod_start."\t".$all_mark[$i]."\t".$all_mark[$i-1]."\n";
#	if($mod_start == 0){next;}
	if($all_mark[$i] ne "*" && $all_mark[$i-1] eq "*"){
		$mod_end = $i-1;
		if($mod_start>0){print O3 "chr\t$mod_start\t$mod_end\tOID_$mod_start-$mod_end\n";}
		$mod_start = $i-1;
	}
	if($all_mark[$i] eq "*" && $all_mark[$i-1] ne "*"){
		$mod_end = $i-1;
		$region_len=$mod_end-$mod_start;
		if($strand eq "-"){$real_pos_use = $real_pos+1;}else{$real_pos_use = $real_pos-$region_len;}
		if($mod_start>0){print O3 "chr\t$mod_start\t$mod_end\tTID_$real_pos_use\n";}
		$mod_start = $i-1;
	}
}
print O2 ">combine\n$final_seq\n";

