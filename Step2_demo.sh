########################
## prepare the following files for Step2 test
# demo/SBDS_ref/combine_mod.fa
# demo/SBDS_ref/combine_mod.info
# demo/SBDS_ref/combine_mod_region.bed
# demo/SBDS_ref/SBDS.bed
# demo/SBDS_ref/SBDSP1.bed
# demo/demo.bam
# demo/demo.bam.bai
########################

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
