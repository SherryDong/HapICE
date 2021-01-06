## prepare the following files for Step2 test
# demo/SBDS_ref/combine_mod.fa
# demo/SBDS_ref/combine_mod.info
# demo/SBDS_ref/combine_mod_region.bed
# demo/SBDS_ref/SBDS.bed
# demo/SBDS_ref/SBDSP1.bed
# demo/demo.bam
# demo/demo.bam.bai

########################
parental_gene=SBDS
pseudo_gene=SBDSP1
total_gene_num=2
result='demo/' ## main directory to save results
combine_ref_dir='demo/SBDS_ref/' ## main directory for combined reference

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
