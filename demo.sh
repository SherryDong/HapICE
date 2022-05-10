
#### Step0: prepare environment parameters
## set parameters
parental_gene=SBDS
pseudo_gene=SBDSP1
total_gene_num=2
## set working directories
current_dir=`pwd`
src_dir=${current_dir}/src
result=${current_dir}/demo ## main directory to save results
combine_ref_dir=${result}/SBDS_ref/ ## main directory for combined reference
db_dir=${current_dir}/db
## user could download from https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta
ref_fa=${db_dir}/Homo_sapiens_assembly19.fasta  ## currently not available, please download 
## user could download from http://www.openbioinformatics.org/annovar/download/hg19_refGene.txt.gz and unzip
ref_annot=${db_dir}/hg19_refGene.txt ## 

#### Step1: prepare gene-specific reference
### Attention: the files under ${db_dir} is currently not available, user need to download files first to process the Step1
## extract fasta files for genes
perl ${src_dir}/gene2fa_combine.pl $parental_gene,$pseudo_gene $combine_ref_dir $ref_fa $ref_annot
### optionals for the first step
## optional I, use GENCODE annotation gtf
## download from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz
## and gunzip into ${db_dir}/gencode.v28lift37.annotation.gtf
# ref_annot_GENCODE='${db_dir}/gencode.v28lift37.annotation.gtf'
# pseudo_gene='SBDSP'
# perl ${src_dir}/gene2fa_combine_GENCODE.pl $true_gene,$pseudo_gene $combine_ref_dir $ref_fa $ref_annot_GENCODE
## optional II, if manually prepare bed file ($ref_annot is not required)
# parental_gene_bed='${result}/SBDS_ref/SBDS.bed'
# pseudo_gene_bed='${result}/SBDS_ref/SBDSP1.bed'
# perl ${src_dir}/gene2fa_combine_manualBed.pl $parental_gene,$pseudo_gene $parental_gene_bed,$pseudo_gene_bed $combine_ref_dir $ref_fa 
###
## perform alignment for parental genes and pseudogenes, generate reference files
muscle -in $combine_ref_dir/combine.fa -out $combine_ref_dir/combine.aln -clw
perl ${src_dir}/common_ref.pl $combine_ref_dir/combine.aln $parental_gene $total_gene_num $combine_ref_dir/$parental_gene.bed $combine_ref_dir/combine_mod.info $combine_ref_dir/combine_mod.fa $combine_ref_dir/combine_mod_region.bed

#### Step2: generate reads-region mapping content matrix for a sample
### Note: we have prepared the following sample files for test: demo/demo.bam,demo/demo.bam.bai; and we have pre-generated reference files under $combine_ref_dir for test
## set input bam files and sample name
bam_file='${result_dir}/demo.bam'
sample_name='demo'
## extract region for parental genes and pseudogenes from bed files
region_parental=`awk -F "\t" {'print "chr"$1":"$2"-"$3'} $combine_ref_dir/$parental_gene.bed | tail -n 1`
region_pseudo=`awk -F "\t" {'print "chr"$1":"$2"-"$3'} $combine_ref_dir/$pseudo_gene.bed | tail -n 1`
## for the target sample, extract aligned reads related with the parental and pseudogenes region
samtools view -h $bam_file $region_parental $region_pseudo >$result/$sample_name.sam
samtools fasta $result/$sample_name.sam >$result/$sample_name.fa
## re-align the reads to the reference files
blat -out=maf $combine_ref_dir/combine_mod.fa $result/$sample_name.fa $result/$sample_name.output.psl
## generate mapping content matrix
perl ${src_dir}/blat_seq2matrix.pl $result/$sample_name.output.psl $combine_ref_dir/combine_mod_region.bed $result/$sample_name.mat

#### Step3: functional/pseudogene haplotype inference and result visualization
### Note: in Step3_demo.R, the test files have been already prepared. 
## read in mapping content matrix and infer the haplotype in the target sample, output detailed mapping statistics with figures
Rscript ${src_dir}/pipeline_draw.R $result/$sample_name.mat $combine_ref_dir/combine_mod_region.bed $combine_ref_dir/combine_mod.info 66459316,66459273,66459256,66459197,66459075,66459073 $result/$sample_name.res.pdf  $result/$sample_name.res.txt


