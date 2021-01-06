#
parental_gene=SBDS
pseudo_gene=SBDSP1
combine_ref_dir='demo/SBDS_ref/' ## main directory for combined reference
## user could download from https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta
ref_fa='db/Homo_sapiens_assembly19.fasta'  ## currently not available, please download 
## user could download from http://www.openbioinformatics.org/annovar/download/hg19_refGene.txt.gz and unzip
ref_annot='db/hg19_refGene.txt' ## currently not available, please download 

# Step1: prepare gene-specific reference
perl src/gene2fa_combine.pl $parental_gene,$pseudo_gene $combine_ref_dir $ref_fa $ref_annot
muscle -in $combine_ref_dir/combine.fa -out $combine_ref_dir/combine.aln -clw
perl src/common_ref.pl $combine_ref_dir/combine.aln $parental_gene $total_gene_num $combine_ref_dir/$parental_gene.bed $combine_ref_dir/combine_mod.info $combine_ref_dir/combine_mod.fa $combine_ref_dir/combine_mod_region.bed
