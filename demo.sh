## prepare the following files for test
# demo/demo.bam
# demo/demo.bam.bai

##
parental_gene=SBDS
pseudo_gene=SBDSP1
total_gene_num=2
result='demo/' ## main directory to save results
combine_ref_dir='demo/SBDS_ref/' ## main directory for combined reference
## user could download from https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta
ref_fa='db/Homo_sapiens_assembly19.fasta'  ## currently not available, please download 
## user could download from http://www.openbioinformatics.org/annovar/download/hg19_refGene.txt.gz and unzip
ref_annot='db/hg19_refGene.txt' ## 

# Step1: prepare gene-specific reference
perl src/gene2fa_combine.pl $parental_gene,$pseudo_gene $combine_ref_dir $ref_fa $ref_annot

### optinals for the first step
## optinal I, use GENCODE annotation gtf
## download from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz
## and gunzip into db/gencode.v28lift37.annotation.gtf
# ref_annot_GENCODE='db/gencode.v28lift37.annotation.gtf'
# pseudo_gene='SBDSP'
# perl src/gene2fa_combine_GENCODE.pl $true_gene,$pseudo_gene $combine_ref_dir $ref_fa $ref_annot_GENCODE
## optional II, if manually prepare bed file ($ref_annot is not required)
# parental_gene_bed='demo/SBDS_ref/SBDS.bed'
# pseudo_gene_bed='demo/SBDS_ref/SBDSP1.bed'
# perl src/gene2fa_combine_manualBed.pl $parental_gene,$pseudo_gene $parental_gene_bed,$pseudo_gene_bed $combine_ref_dir $ref_fa 
###

muscle -in $combine_ref_dir/combine.fa -out $combine_ref_dir/combine.aln -clw
perl src/common_ref.pl $combine_ref_dir/combine.aln $parental_gene $total_gene_num $combine_ref_dir/$parental_gene.bed $combine_ref_dir/combine_mod.info $combine_ref_dir/combine_mod.fa $combine_ref_dir/combine_mod_region.bed

# Step2: generate reads-region mapping content matrix for a sample
bam_file='demo/demo.bam'
sample_name='demo'

region_parental=`awk -F "\t" {'print "chr"$1":"$2"-"$3'} $combine_ref_dir/$parental_gene.bed | tail -n 1`
region_pseudo=`awk -F "\t" {'print "chr"$1":"$2"-"$3'} $combine_ref_dir/$pseudo_gene.bed | tail -n 1`
# for each sample
samtools view -h $bam_file $region_parental $region_pseudo >$result/$sample_name.sam
samtools fasta $result/$sample_name.sam >$result/$sample_name.fa
blat -out=maf $combine_ref_dir/combine_mod.fa $result/$sample_name.fa $result/$sample_name.output.psl
perl src/blat_seq2matrix.pl $result/$sample_name.output.psl $combine_ref_dir/combine_mod_region.bed $result/$sample_name.mat

# Step3: functional/pseudogene haplotype inference and result visualization
Rscript src/pipeline_draw.R $result/$sample_name.mat $combine_ref_dir/combine_mod_region.bed $combine_ref_dir/combine_mod.info 66459316,66459273,66459256,66459197,66459075,66459073 $result/$sample_name.res.pdf  $result/$sample_name.res.txt


