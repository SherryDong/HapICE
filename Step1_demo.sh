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
