import os
import pandas as pd
import numpy as np

id = ['CNE1', 'KM4', 'M368']
depth = pd.read_table('/project/jnovembre/chichun/ancient-imputation-experiment/data/input/sample_depth.tbl', sep=' ', header=None)

## create a dataframe of sample fraction in each sample
## a tiny bit convoluted and to be simplified
coverage = [0.1 , 0.5 , 1 , 2]
dp = {}
for s in id:
    dp[s] = list(coverage/depth[depth.iloc[:,0]==s].iloc[:,1].values[0])
dp = pd.DataFrame(dp)  
coverage = [str(c) for c in coverage]
dp['coverage']=coverage


## start of the rules
rule ref_all:
    input:
        expand("data/bam/1240k/downsampled/{id}_allmerged_rmdup_reheader_1240k_{coverage}x.trimmed.bam", id = id , coverage = coverage)
        

rule down_sample_bams:
    '''
    downsample bams from a fraction of reads (default ~5x, 2x, 1x, 0.5x, 0.1x)
    '''
    input:
        "data/bam/1240k/{id}_allmerged_rmdup_reheader_1240k.trimmed.bam"
    output:
        "data/bam/1240k/downsampled/{id}_allmerged_rmdup_reheader_1240k_{coverage}x.trimmed.bam"
    threads: 12
    run:
        # get the read sampling fraction from the dataframe with sample id and target coverage
        frac = float(dp[dp['coverage']== wildcards.coverage ][wildcards.id].values[0])
        shell("samtools view -bs {frac} {input} > {output}")
        shell("samtools index {output}")
          
 
