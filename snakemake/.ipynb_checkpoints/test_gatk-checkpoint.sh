#!/bin/sh
module load java

bam='/project/jnovembre/chichun/ancient-imputation-experiment/data/input/bam.wgs.list'
ref_1kg='/project/jnovembre/chichun/ancient-imputation-experiment/data/input/allele/wgs_allele.chr2.vcf.gz'
dbsnp='/project/jnovembre/chichun/ancient-imputation-experiment/data/input/dbsnp_138.b37.vcf'
ref_genome='/project/jnovembre/chichun/ancient-imputation-experiment/data/input/hs37d5.fa'

java -Xmx30g -jar utils/GenomeAnalysisTK.jar \
      -T UnifiedGenotyper \
      -I ${bam} \
      -L ${ref_1kg} \
      -R ${ref_genome} \
      --genotype_likelihoods_model SNP \
      -mbq 30 \
      -stand_call_conf 0.0 \
      --genotyping_mode GENOTYPE_GIVEN_ALLELES \
      -alleles ${ref_1kg} \
      --allSitePLs \  
      --dbsnp ${dbsnp} \
      -rf BadCigar \
      --output_mode EMIT_ALL_SITES \
      -o gatk.wgs.chr2.vcf.gz
