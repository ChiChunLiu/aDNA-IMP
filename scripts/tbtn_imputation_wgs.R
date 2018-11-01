.libPaths(c("/home/chichun/R_libs","/home/chichun/R/x86_64-pc-linux-gnu-library/3.4"))
library(GeneImp)

setwd('/project/jnovembre/chichun/ancient-imputation-experiment')
sample.dir = 'data/vcf/wgs'
reference.dir = 'data/1kg_reference/min_mac5'

args = commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])
j = as.numeric(args[2])
k = as.numeric(args[3])

options('vcf.yieldSize'=1E5)

imputevcf(vcfname = paste0(sample.dir,"/gatk.wgs.masked.merged.chr", i, ".vcf.gz"),
          ref.vcfname = paste0(reference.dir,"/1KG_HA.biallelic.mac5.chr", i, ".vcf.gz"),
          klthresh = j,
	        flanksize = 0.5,
          filtermethod = "pairrand",
	        numfilterhaps = 200,
          maxjobs = k,
	        write.dir = paste0(sample.dir,"/imputed_data_kl", j),
          temp.dir = NULL,
	        verbose = 0,
	        diagnostics = FALSE)

