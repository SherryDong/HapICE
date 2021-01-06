#!/usr/bin/perl -w

$psl = $ARGV[0];
$mod_region = $ARGV[1];
###############
#chr 98  119 OID_98-119
open M,$mod_region or die $!;
while(<M>){
	chomp;
	@a = split "\t";
	if($a[3] =~ /OID/){next;}
	$region_start{$a[3]} = $a[1];
	$region_end{$a[3]} = $a[2];
	foreach $s ($a[1]..$a[2]){
		$pos2region{$s} = $a[3];
	}
}
###############
open O,$psl or die $!;
## s combine                                 7526 150 + 8107 gttaacagtggacacagcaaggctttcc...
## s A01195:14:HHYYGDSXY:2:1514:5710:21778/1    0 150 -  150 gttaacagtggacacagcaaggctttc...
while(<O>){
	chomp;
	$line = $_;
	if(/score=/){next;}
	if($line eq ""){next;}
	if(/^#/){next;}
	if(/^s combine\s+(\d+)\s+(\d+)\s+([+|-])\s+(\d+)\s+(.*)/){
		$ref_start = $1;
		$ref_len = $2;
		#$ref_strand = $3;
		$ref_end = $ref_start+$ref_len-1;
		#$ref_seq = $5;
	}else{
		if(/^s (.*)\s+(\d+)\s+(\d+)\s+([+|-])\s+(\d+)\s+(.*)/){
			$name = $1;
			#$start = $2;
			#$len = $3;
			#$strand = $4;
			#$end = $5;
			$seq = $6;
		}
		##
		$seq_name = $name;
		if($name =~/(.*)\/[1|2]/){
			$seq_name = $1;
		}
		## get TID
		undef(%each_region);
		foreach $s ($ref_start..$ref_end){
			if($pos2region{$s}){
				$region = $pos2region{$s};
				$each_region{$region} = 1;
			}
		}
		foreach $region (keys %each_region){
				$rs = $region_start{$region};
				$re = $region_end{$region};
				$len = $re-$rs;
				$ss = $rs-$ref_start+1;
				$target_seq = substr($seq,$ss,$len); 
				#print "$region\t$rs\t$re\t$ref_start\t$ref_end\t$ss\t$len\t$target_seq\t$.\n";
				if($target_seq eq ""){next;}
				$final{$seq_name}{$region} = uc($target_seq);
				$all_region{$region} = 1;
		}
	}
}
close O;

#### output
open O,">$ARGV[2]";
@all_region = sort(keys %all_region);
$n = join "\t",@all_region;
print O "seq_name\t$n\n";
foreach $seq_name (keys %final){
	undef(@n);
	foreach $region (@all_region){
		if($final{$seq_name}{$region}){
			push(@n,$final{$seq_name}{$region});		
		}else{
			push(@n,".");
		}
	}
	$n = join "\t",@n;
	print O "$seq_name\t$n\n";
}
close O;
