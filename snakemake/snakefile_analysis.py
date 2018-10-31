configfile: '/project/jnovembre/chichun/ancient-imputation-experiment/config/config.yaml'
shell.prefix("module load java; ")

chr = [str(c) for c in list(range(1,23))]

rule all:
    input:
        expand('{}/{}'.format(config['path']['vcf_output'], 'f50k/gatk.f50k.chr{CHR}.vcf.gz'), CHR = chr),
        expand('{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.f50k.chr{CHR}.vcf.gz'), CHR=chr),
        expand('{}/{}'.format(config['path']['vcf_output'], 'f50k/KS20_KS25.S143.f50k.chr{CHR}.gpaveraged.vcf.gz'), CHR=chr)
        
rule subset_f50k_empty_vcf:
    '''
    subset functional 50k reference vcf
    '''
    input:
        vcf = '{}/{}'.format(config['path']['data_input'], 'allele/wgs_allele.chr{CHR}.vcf.gz'),
        pos = '{}/{}'.format(config['path']['data_input'], 'functional50k.chr{CHR}.pos')
    output:
        '{}/{}'.format(config['path']['data_input'], 'allele/f50k_allele.chr{CHR}.vcf.gz')
    run:
        shell('bcftools view -R {input.pos} {input.vcf} -Oz -o {output} && tabix -p vcf {output}')
        
rule subset_f50k_reference_vcf:
    '''
    subset functional 50k reference vcf
    '''
    input:
        vcf = '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.chr{CHR}.vcf.gz'),
        pos = '{}/{}'.format(config['path']['data_input'], 'functional50k.chr{CHR}.pos')
    output:
        '{}/{}'.format(config['path']['ref_output'], '1KG_HA.biallelic.mac5.f50k.chr{CHR}.vcf.gz')
    run:
        shell('bcftools view -R {input.pos} {input.vcf} -Oz -o {output} && tabix -p vcf {output}')

        
rule f50k_capture_call:
    '''
    call genotypes KS20_KS25, S143_S173
    '''
    input:
        bam = '{}/{}'.format(config['path']['data_input'], 'bam.KS20_KS25.S143_S173.list'),
        ref_1kg = '{}/{}'.format(config['path']['data_input'], 'allele/f50k_allele.chr{CHR}.vcf.gz'),
        dbsnp = '{}/{}'.format(config['path']['ref_genome'], 'dbsnp_138.b37.vcf'),
        ref_genome = '{}/{}'.format(config['path']['ref_genome'], 'hs37d5.fa')
    output:
        '{}/{}'.format(config['path']['vcf_output'], 'f50k/gatk.f50k.chr{CHR}.vcf.gz')
    run:
        shell('java -Xmx2g -jar utils/GenomeAnalysisTK.jar ' \
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
        
rule average_f50k_ancient:
    input:
        vcf1 = '{}/{}'.format(config['path']['vcf_output'], 'wgs/imputed_data_kl15/all.imputed.gatk.wgs.masked.merged.KS20_KS25.S143.chr{CHR}.kl15.flank50.haps200.pairrand.vcf.gz'),
        vcf2 = '{}/{}'.format(config['path']['vcf_output'], 'wgs/imputed_data_kl20/all.imputed.gatk.wgs.masked.merged.KS20_KS25.S143.chr{CHR}.kl20.flank50.haps200.pairrand.vcf.gz'),
        vcf3 = '{}/{}'.format(config['path']['vcf_output'], 'wgs/imputed_data_kl25/all.imputed.gatk.wgs.masked.merged.KS20_KS25.S143.chr{CHR}.kl25.flank50.haps200.pairrand.vcf.gz')
    output:
        wgs = '{}/{}'.format(config['path']['vcf_output'], 'wgs/KS20_KS25.S143.chr{CHR}.gpaveraged.vcf.gz')
    run:
        shell('python scripts/genoAverage.py -a {input.vcf1} -b {input.vcf2} -c {input.vcf3} -o {output.wgs} && tabix -p vcf {output.wgs}')
        

rule filter_f50k_snps:
    input:
        pos = '{}/{}'.format(config['path']['data_input'], 'functional50k.chr{CHR}.pos'),
        wgs = '{}/{}'.format(config['path']['vcf_output'], 'wgs/KS20_KS25.S143.chr{CHR}.gpaveraged.vcf.gz')
    output:
         '{}/{}'.format(config['path']['vcf_output'], 'f50k/KS20_KS25.S143.f50k.chr{CHR}.gpaveraged.vcf.gz')
    run:
        shell('bcftools view -R {input.pos} {input.wgs} -Oz -o {output} && tabix -p vcf {output}')
