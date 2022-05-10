# HapICE V1.0

HapICE (Haplotype Inference for Pseudogene-mediated conversion events based on short-read next-generation sequencing data) is a solution started from the general short-read NGS BAM files and dedicated to inferring the genomic haplotypes of functional/pseudogene pairs that have highly-homologous sequences and gene conversion events frequently occurred, along with result visualization including inferred haplotype, the proportion of gene recombination events between two neighboring informative mutations, and specific reads information. 


# Require softwares

```seqtk >= 1.3-r106```

```bedtools >= v2.27.1```

```muscle >= v3.8.31```

```samtools >= v1.8```

```blat >= v.36```

```R >= 3.6.0```

```perl >= 5```

# Require files

In Step1, users need to download database files under the instruction and put them under db/ directory (detailed steps check demo.sh or Step1_demo.sh). This step is not required to test Step2_demo.sh and Step3_demo.R. In Step2_demo.sh and Step3_demo.R, the test files have been prepared.  

# Demo

HapICE performs haplotype inference followed by three main steps: (1) prepare gene-specific reference, (2) generate reads-region mapping content matrix, and (3) functional/pseudogene haplotype inference and result visualization 

Suggest users to test files under the main working directory, do not change directory to src/ or result/.

### Demos script for all of the process.
demo.sh

### Demos script for each of the three steps (each script could run separately). 

+ Part I: prepare gene-specific reference

  Step1_demo.sh

+ Part II: generate reads-region mapping content matrix

  Step2_demo.sh
  
  # 

+ Part III: functional/pseudogene haplotype inference and result visualization 

  Step3_demo.R
  
  Rscript src/pipeline_draw.R

# docker image
docker pull jargene/hapice:1.0

docker run -d -t -i docker.io/jargene/hapice:1.0 "/bin/bash"

# Reference
medRxiv: https://medrxiv.org/cgi/content/short/2021.08.22.21262444v1
