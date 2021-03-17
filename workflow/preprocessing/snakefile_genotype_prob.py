configfile: '/project/jnovembre/chichun/ancient-imputation-experiment/config/config.yaml'
shell.prefix("module load java; ")
# create bam lists as inputs for gatk unified genotyper.
#os.system('ls data/bam/1240k/downsampled/*.bam > data/input/1240k_downsampled_bam.list')
#os.system('ls data/bam/1240k/*.bam >> data/input/1240k_downsampled_bam.list')

#os.chdir('/project/jnovembre/chichun/ancient-imputation-experiment')
#os.system('ls /project/jnovembre/chichun/ancient-tbtn-analysis/data/bam/wgs_read/*trimmed.reheader.bam > data/input/bam.wgs.list')
#with open("data/input/bam.wgs.list", "a") as myfile:
#    myfile.write('{}/{}'.format(config['path']['bam'], 'sequence/LBK.hg19_1000g.2x.bam\n'))
#    myfile.write('{}/{}'.format(config['path']['bam'], 'sequence/LP6005443-DNA_E09.2x.srt.aln.bam\n'))
#    myfile.write('{}/{}'.format(config['path']['bam'], 'sequence/LP6005442-DNA_G01.2x.srt.aln.bam\n'))
    
#os.system('ls /project/jnovembre/chichun/ancient-imputation-experiment/data/bam/1240k/*.capture.reheader.trimmed.bam > data/input/bam.capture.list')

chr = [str(c) for c in list(range(1,23))]

rule all:
    input:
        #expand('{}/{}'.format(config['path']['vcf_output'], '1240k/gatk.1240k.chr{CHR}.vcf.gz'),CHR=chr),
        expand('{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.masked.merged.chr{CHR}.vcf.gz'), CHR=chr)
        #expand('{}/{}'.format(config['path']['vcf_output'], '1240k/mpileup.1240k.chr{CHR}.vcf.gz'), CHR=chr),
        #expand('{}/{}'.format(config['path']['vcf_output'], 'wgs/mpileup.wgs.chr{CHR}.vcf.gz'), CHR=chr),
        #expand('{}/{}'.format(config['path']['vcf_output'], '1240k/rdraw.1240k.chr{CHR}.vcf.gz'), CHR=chr),
        #expand('{}/{}'.format(config['path']['vcf_output'], 'wgs/rdraw.wgs.chr{CHR}.vcf.gz'), CHR=chr),
        #expand('{}/{}'.format(config['path']['vcf_output'], '1240k/aaca.trimmed.1240k.randomdraw.chr{CHR}.vcf.gz'), CHR=chr)        

rule comp_lik_gatk_wgs:
    '''
    compute genotype likelihood from gatk, wgs data
    '''
    input:
        # samples: 18 wgs aACA + 1 LBK + 2 SGDP
        bam = '{}/{}'.format(config['path']['data_input'], 'bam.wgs.list'),
        ref_1kg = '{}/{}'.format(config['path']['data_input'], 'allele/wgs_allele.chr{CHR}.vcf.gz'),
        dbsnp = '{}/{}'.format(config['path']['ref_genome'], 'dbsnp_138.b37.vcf'),
        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
    output:
        '{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.chr{CHR}.vcf.gz')
    run:
        shell('java -Xmx25g -jar utils/GenomeAnalysisTK.jar ' \
              '-T UnifiedGenotyper ' \
              '-I {input.bam} ' \
              '-L {input.ref_1kg} ' \
              '-R {input.ref_genome} ' \
              '--genotype_likelihoods_model SNP ' \
              '-mbq 30 ' \
              '-stand_call_conf 0.0 ' \
              '--genotyping_mode GENOTYPE_GIVEN_ALLELES ' \
              '-alleles {input.ref_1kg} ' \
              '--allSitePLs ' \  
              '--dbsnp {input.dbsnp} ' \
              '-rf BadCigar ' \
              '--output_mode EMIT_ALL_SITES ' \
              '-o {output}')
                                 
rule comp_lik_gatk_1240k:
    '''
    compute genotype likelihood from gatk, 1240k sites only
    '''
    input:
        # samples: 23 aACA with capture
        bam = '{}/{}'.format(config['path']['data_input'], 'bam.capture.list'),
        ref_1kg = '{}/{}'.format(config['path']['data_input'], 'allele/1240k_allele.chr{CHR}.vcf.gz'),
        dbsnp = '{}/{}'.format(config['path']['ref_genome'], 'dbsnp_138.b37.vcf'),
        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
    output:
        '{}/{}'.format(config['path']['vcf_output'], '1240k/gatk.1240k.chr{CHR}.vcf.gz')        
    run:
        shell('java -Xmx10g -jar utils/GenomeAnalysisTK.jar ' \
              '-T UnifiedGenotyper ' \
              '-I {input.bam} ' \
              '-L {input.ref_1kg} ' \
              '-R {input.ref_genome} ' \
              '--genotype_likelihoods_model SNP ' \
              '--allSitePLs ' \
              '-mbq 30 ' \
              '-stand_call_conf 0.0 ' \
              '--genotyping_mode GENOTYPE_GIVEN_ALLELES ' \
              '-alleles {input.ref_1kg} ' \
              '--dbsnp {input.dbsnp} '\
              '-rf BadCigar ' \
              '--output_mode EMIT_ALL_SITES ' \
              '-o {output} ')

rule pileup_wgs:
    '''
    bcftools mpileup for drawing pseudo haploid genome
    '''
    input:
        bam = '{}/{}'.format(config['path']['data_input'], 'bam.wgs.list'),
        tsv = '{}/{}'.format(config['path']['data_input'], 'pos/wgs.chr{CHR}.tsv.gz'),
        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
    output:
        '{}/{}'.format(config['path']['vcf_output'], 'wgs/mpileup.wgs.chr{CHR}.vcf.gz')
    run:
        shell('bcftools mpileup --ignore-RG -B -q30 -Q30 ' \
              '-T {input.tsv} '\
              '-f {input.reference_genome} ' \
              '-b {input.bam} ' \
              '-a AD,DP | ' \
              'bcftools norm -m-both -Oz -o {output}' )
        shell('tabix -p vcf {output}')
              
rule pileup_1240k:
    '''
    bcftools mpileup for drawing pseudo haploid genome
    '''
    input:
        bam = '{}/{}'.format(config['path']['data_input'], 'bam.capture.list'),
        tsv = '{}/{}'.format(config['path']['data_input'], 'pos/1240k2ref.chr{CHR}.tsv.gz'),
        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
    output:
        '{}/{}'.format(config['path']['vcf_output'], '1240k/mpileup.1240k.chr{CHR}.vcf.gz')  
    run:
        shell('bcftools mpileup --ignore-RG -B -q30 -Q30 ' \
              '-T {input.tsv} '\
              '-f {input.reference_genome} ' \
              '-b {input.bam} ' \
              '-a AD,DP | ' \
              'bcftools norm -m-both -Oz -o {output}' )
        shell('tabix -p vcf {output}')

rule random_pseudo_haploid:
    '''
    Draw random allele with my messy but fast python script
    with bcftools mpileup file as input
    '''
    input:
        snp_1240k = '{}/{}'.format(config['path']['data_input'], 'pos/1240k2ref.chr{CHR}.snp'),
        snp_wgs = '{}/{}'.format(config['path']['data_input'], 'pos/wgs.chr{CHR}.snp'),
        sample_1240k = '{}/{}'.format(config['path']['data_input'], 'sample.1240k.txt'),
        sample_wgs = '{}/{}'.format(config['path']['data_input'], 'sample.wgs.txt'),
        pileup_1240k = '{}/{}'.format(config['path']['vcf_output'], '1240k/mpileup.1240k.chr{CHR}.vcf.gz'),
        pileup_wgs = '{}/{}'.format(config['path']['vcf_output'], 'wgs/mpileup.wgs.chr{CHR}.vcf.gz')
    output:
        vcf_1240k = '{}/{}'.format(config['path']['vcf_output'], '1240k/rdraw.1240k.chr{CHR}.vcf.gz'),
        vcf_wgs = '{}/{}'.format(config['path']['vcf_output'], 'wgs/rdraw.wgs.chr{CHR}.vcf.gz')
    run:
        shell('scripts/mpileup2vcf.py -s {input.snp_1240k} -i {input.sample_1240k} -p {input.pileup_1240k} -o {output.vcf_1240k}')
        shell('scripts/mpileup2vcf.py -s {input.snp_1240k} -i {input.sample_wgs} -p {input.pileup_wgs} -o {output.vcf_wgs}')

rule mask_transition_PL:
    '''
    masking pred-scale likelihood to be 0,0,0 to reduce bias from damage
    '''
    input:
        vcf_wgs = '{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.chr{CHR}.vcf.gz')
    output:
        vcf_wgs = '{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.masked.chr{CHR}.vcf.gz')
    run:
        shell('python scripts/mask_transition_PL.py -v {input.vcf_wgs}  -o {output.vcf_wgs}')
        shell('tabix -p vcf {output.vcf_wgs}')

rule rename_unmasked_sgdpLBK:
    '''
    Imputation is performed on both masked and unmasked Naxi, Yi, LBK samples, 
    so we need to create different sample names
    '''
    input:
        vcf = '{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.chr{CHR}.vcf.gz'),
        sample = '{}/{}'.format(config['path']['data_input'], 'sgdp_LBK_sample.txt')
    output:
        '{}/{}'.format(config['path']['vcf_output'], 'wgs/tmp.chr{CHR}.vcf.gz') 
    run:
        shell('bcftools view -S {input.sample} {input.vcf} -Oz -o {output}')
        shell('scripts/rename_vcf_sample.sh {output} _unmasked')

rule merge_renamedSGDPLBK_masked:
    '''
    Merged renamed Naxi, Yi, LBK with all masked samples
    '''
    input:
        vcf_unmasked = '{}/{}'.format(config['path']['vcf_output'], 'wgs/tmp.chr{CHR}.vcf.gz'),
        vcf_masked = '{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.masked.chr{CHR}.vcf.gz')
    output:
        '{}/{}'.format(config['path']['vcf_output'], 'wgs/gatk.wgs.masked.merged.chr{CHR}.vcf.gz')       
    run:
        shell('bcftools merge {input.vcf_unmasked} {input.vcf_masked} -Oz -o {output}')
        shell('tabix -p vcf {output}')
        shell('rm {input.vcf_unmasked}*')
              
#rule comp_lik_gatk_downsampled:
#    '''
#    compute genotype likelihood from gatk, all downsampled data
#    '''
#    input:
#        # samples: downsampled aACA 
#        bam = 
#        ref_1kg = '{}/{}'.format(config['path']['ref_input'],'1KG_HA.biallelic.mac5.chr{CHR}.vcf.gz',
#        dbsnp = '{}/{}'.format(config['path']['ref_genome'], 'dbsnp_138.b37.vcf')
#        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
#    output:
#        '{}/{}'.format(config['path']['vcf_output'], 'wgs/aaca_trimmed_downsampled.chr{CHR}.vcf.gz')
#    run:
#        shell('java -Xmx10g -jar utils/GenomeAnalysisTK.jar ' \
#              '-T UnifiedGenotyper ' \
#              '-I {input.bam} ' \
#              '-L {input.ref_1kg} '\
#              '-R {input.ref_genome} ' \
#              '--genotype_likelihoods_model SNP ' \
#              '-mbq 30 ' \
#              '-stand_call_conf 0.0 ' \
#              '--genotyping_mode GENOTYPE_GIVEN_ALLELES ' \
#              '-alleles {input.ref_1kg} ' \
#              '--allSitePLs ' \  
#              '--dbsnp {input.dbsnp} ' \
#              '-rf BadCigar ' \
#              '--output_mode EMIT_ALL_SITES ' \
#              '-o {output}')

#rule random_pseudo_hap_mathieson:
#    '''
#    randomly draw psdeudo haploid data with Mathieson gdc python script
#    see gdc repo at: https://github.com/mathii
#    ''' 
#    input: 
#        bam_list = 'data/input/downsample_randomdraw_bam.list',
#        snp_list = 'data/input/downsample_randomdraw_chr{CHR}.snp'
#    output:
#        'data/vcf/1240k/pseudo_haploid/aaca_1240k_chr{CHR}_randomdraw.vcf'
#    run:
#        shell('python utils/gdc/apulldown.py ' \
#              '-b {input.bam_list} -s {input.snp_list} -q 30 -i 30 ' \
#              '-o vcf > {output}')
