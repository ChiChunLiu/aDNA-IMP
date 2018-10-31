import os
import pandas as pd
import numpy as np
configfile: '/project/jnovembre/chichun/ancient-imputation-experiment/config/config.yaml'
    
# samples for downsample experiment
# these three samples are sequenced and captured to intermedium coverage 
id = ['CNE1', 'KM4', 'M368']
depth = pd.read_table('{}/{}'.format(config['path']['data_input'],'sample_depth.tbl'),
                      sep=' ', header=None)

# samples for 1240k imputation
df_1240k = pd.read_table('{}/{}'.format(config['path']['data_input'],'capture_1240k.sample.list'),
                         sep = '\t', header = None, names = ['id', 'rg'])
id_1240k = list(df_1240k['id'].values)

id_wgs = os.listdir('{}/{}'.format(config['path']['merged_bam'], 'wgs_read'))
id_wgs = [id.split('.trimmed.bam')[0] for id in id_wgs if id.endswith('.trimmed.bam')]

## create a dataframe of sample fraction in each sample
## a bit convoluted and to be simplified
coverage = [0.1 , 0.5 , 1 , 2]
dp = {}
for s in id:
    dp[s] = list(coverage/depth[depth.iloc[:,0]==s].iloc[:,1].values[0])
dp = pd.DataFrame(dp)  
coverage = [str(c) for c in coverage]
dp['coverage']=coverage


# rules
rule all:
    input:
        '{}/{}'.format(config['path']['bam'], 'sequence/LBK.hg19_1000g.2x.bam'),
        '{}/{}'.format(config['path']['bam'], 'sequence/LP6005443-DNA_E09.2x.srt.aln.bam'),
        '{}/{}'.format(config['path']['bam'], 'sequence/LP6005442-DNA_G01.2x.srt.aln.bam'),
        expand('{}/{}'.format(config['path']['merged_bam'], 'wgs_read/{id_wgs}.trimmed.reheader.bam'), id_wgs = id_wgs),
        expand('{}/{}'.format(config['path']['bam'], '1240k/{id_1240k}.capture.reheader.trimmed.bam'), id_1240k = id_1240k)
        #expand('{}/{}'.format(config['path']['bam'], '1240k/downsampled/{id}.capture.{coverage}x.trimmed.bam'), id = id , coverage = coverage)


rule downsample_bam_LBK:
    '''
    downsample LBK to 2x 
    '''
    input:
	    '{}/{}'.format(config['path']['bam'], 'sequence/LBK.hg19_1000g.bam')
    output:
	    '{}/{}'.format(config['path']['bam'], 'sequence/LBK.hg19_1000g.2x.bam')
    run:
        shell('samtools view -bs 0.1047 {input} > {output}')
        shell('samtools index {output}')
        
rule downsample_bam_Naxi:
    '''
    downsample Naxi to 2x 
    '''
    input:
        '{}/{}'.format(config['path']['sgdp_bam'], 'LP6005443-DNA_E09.srt.aln.bam')
    output:
        '{}/{}'.format(config['path']['bam'], 'sequence/LP6005443-DNA_E09.2x.srt.aln.bam')
    run:
        shell('samtools view -bs 0.0484 {input} > {output}')
        shell('samtools index {output}')
        
rule downsample_bam_Yi:
    '''
    downsample Yi to 2x
    '''
    input:
        '{}/{}'.format(config['path']['sgdp_bam'], 'LP6005442-DNA_G01.srt.aln.bam')
    output:
        '{}/{}'.format(config['path']['bam'], 'sequence/LP6005442-DNA_G01.2x.srt.aln.bam')
    run:
        shell('samtools view -bs 0.048 {input} > {output}')
        shell('samtools index {output}')
                
rule filter_1240k_read:
    '''
    filter 1240K capture reads from read group information
    '''
    input:
        bam = '{}/{}'.format(config['path']['merged_bam'], '{id_1240k}_allmerged_rmdup_reheader.trimmed.bam'),
        rg = '{}/{}'.format(config['path']['data_input'], 'readgroup/{id_1240k}.rg')
    output:
	    '{}/{}'.format(config['path']['bam'], '1240k/{id_1240k}.capture.trimmed.bam')
    run:
        shell('samtools view -hbR {input.rg} {input.bam} > {output}')
        shell("samtools index {output}")

rule downsample_bam_1240k:
    '''
    downsample bams from a fraction of reads (default ~5x, 2x, 1x, 0.5x, 0.1x)
    '''
    input:
        '{}/{}'.format(config['path']['bam'], '1240k/{id}.capture.trimmed.bam')
    output:
        '{}/{}'.format(config['path']['bam'], '1240k/downsampled/{id}.capture.{coverage}x.trimmed.bam')
    run:
        # get the read sampling fraction from the dataframe with sample id and target coverage
        frac = float(dp[dp['coverage'] == wildcards.coverage][wildcards.id].values[0])
        shell("samtools view -bs {frac} {input} > {output}")
        shell("samtools index {output}")
        
rule rename_bam_sm_wgs:
    '''
    rename sm field in bam files
    '''
    input:
        '{}/{}'.format(config['path']['merged_bam'], 'wgs_read/{id_wgs}.trimmed.bam')
    output:
	    '{}/{}'.format(config['path']['merged_bam'], 'wgs_read/{id_wgs}.trimmed.reheader.bam')
    run:
        sample = wildcards.id_wgs
        shell('bash /project/jnovembre/chichun/ancient-imputation-experiment/utils/rename_header.sh ' \
              '{input} {sample} {output}')
        
rule rename_bam_sm_1240k:
    '''
    rename sm field in bam files
    '''
    input:
        '{}/{}'.format(config['path']['bam'], '1240k/{id_1240k}.capture.trimmed.bam')
    output:
	    '{}/{}'.format(config['path']['bam'], '1240k/{id_1240k}.capture.reheader.trimmed.bam')
    run:
        sample = wildcards.id_1240k
        shell('bash /project/jnovembre/chichun/ancient-imputation-experiment/utils/rename_header.sh ' \
              '{input} {sample} {output}')