configfile: '/project/jnovembre/chichun/ancient-imputation-experiment/config/config.yaml'
    
chr = [str(c) for c in list(range(1,23))]

rule all:
    input:
        expand('{}/{}'.format(config['path']['data_input'], 'pos/1240k2ref.chr{CHR}.tsv.gz'), CHR = chr),
        expand('{}/{}'.format(config['path']['data_input'], 'pos/wgs.chr{CHR}.tsv.gz'), CHR = chr),
        expand('{}/{}'.format(config['path']['data_input'], 'pos/1240k2ref.chr{CHR}.snp'), CHR = chr),
        expand('{}/{}'.format(config['path']['data_input'], 'pos/wgs.chr{CHR}.snp'), CHR = chr),
        expand('{}/{}'.format(config['path']['data_input'], 'allele/1240k_allele.chr{CHR}.vcf.gz'), CHR = chr),
        expand('{}/{}'.format(config['path']['data_input'], 'allele/wgs_allele.chr{CHR}.vcf.gz'), CHR = chr)

rule norm_biallelic_mac5:
    '''
    merge multiallelic records and keep biallelic SNPs with minimal allele count 5
    '''
    input:
        vcf = '{}/{}'.format(config['path']['ref_input'], '1KG_HA_chr{CHR}.vcf.gz'),
        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
    output:
        '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.chr{CHR}.vcf.gz')
    run:
        shell('bcftools norm -m+both {input.vcf} | bcftools norm -f {input.ref_genome} | ' \
              'bcftools view -c5:minor -m2 -M2 -v snps -Oz -o {output}')
        shell('tabix -p vcf {output}')


rule filter_1240k_sites:
    '''
    filter 1240k sites in reference panel
    '''
    input:
        vcf = '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.chr{CHR}.vcf.gz'),
        pos = '{}/{}'.format(config['path']['data_input'], '1240k.chr{CHR}.pos')
    output:
        '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.1240k.chr{CHR}.vcf.gz')
    run:
        shell('bcftools view -R {input.pos} {input.vcf} -Oz -o {output}')
        shell('tabix -p vcf {output}')
        
rule drop_genotype:
    '''
    get 1240k position present in reference panel. Also creates vcf 
    position files for both wgs and 1240k sites in reference panel,
    dropping individual GT info.
    
    Note that, querying position may result in duplicate records. Here 
    we're fine since they are all biallelic records.
    '''
    input:
        vcf_1240k = '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.1240k.chr{CHR}.vcf.gz'),
        vcf_wgs = '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.chr{CHR}.vcf.gz')
    output:
        vcf_1240k = '{}/{}'.format(config['path']['data_input'], 'allele/1240k_allele.chr{CHR}.vcf.gz'),
        vcf_wgs = '{}/{}'.format(config['path']['data_input'], 'allele/wgs_allele.chr{CHR}.vcf.gz')
    run:
        # -G/--drop-genotypes is used as input in GATK unified genotyper.
        shell('bcftools view -G {input.vcf_1240k} -Oz -o {output.vcf_1240k}')
        shell('bcftools view -G {input.vcf_wgs} -Oz -o {output.vcf_wgs}')
        shell('tabix -p vcf {output.vcf_1240k}')
        shell('tabix -p vcf {output.vcf_wgs}')
        
rule final_tsv_snp:
    '''
    convert reference vcf into snp and tsv file
    '''
    input:
        vcf_1240k = '{}/{}'.format(config['path']['data_input'], 'allele/1240k_allele.chr{CHR}.vcf.gz'),
        vcf_wgs = '{}/{}'.format(config['path']['data_input'], 'allele/wgs_allele.chr{CHR}.vcf.gz')
    output:
        tsv_1240k = '{}/{}'.format(config['path']['data_input'], 'pos/1240k2ref.chr{CHR}.tsv.gz'),
        tsv_wgs = '{}/{}'.format(config['path']['data_input'], 'pos/wgs.chr{CHR}.tsv.gz'),
        snp_1240k = '{}/{}'.format(config['path']['data_input'], 'pos/1240k2ref.chr{CHR}.snp'),
        snp_wgs = '{}/{}'.format(config['path']['data_input'], 'pos/wgs.chr{CHR}.snp')
    run:
        shell("bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' {input.vcf_1240k} | " \
              "bgzip -c > {output.tsv_1240k} && tabix -s1 -b2 -e2 {output.tsv_1240k}")
        shell("bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' {input.vcf_wgs} | " \
              "bgzip -c > {output.tsv_wgs} && tabix -s1 -b2 -e2 {output.tsv_wgs}")
        shell("bcftools query -f'%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\%REF\\%ALT\\n' " \
              "{input.vcf_1240k} > {output.snp_1240k}")
        shell("bcftools query -f'%CHROM:%POS:%REF:%ALT\\t%CHROM\\t%POS\\%REF\\%ALT\\n' " \
              "{input.vcf_wgs} > {output.snp_wgs}")
        