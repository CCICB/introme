process variant_info {
    container "${params.variant_info_docker_container}"
    beforeScript 'echo Starting variant_info'
    afterScript  'echo Completed variant_info'
    publishDir (path: "${params.outdir}/variant_info")

    input:
        path input_vcf
        path gtf

    output:
        path "${params.prefix}.variant_info.filtered.vcf.gz", emit: variant_info
        path "${params.prefix}.variant_info.filtered.vcf.gz.tbi", emit: variant_info_tbi
        path "${params.prefix}.variant_info.filtered.rmanno.vcf.gz", emit: variant_info_rmanno

    script:
    """
    gunzip -c $gtf | sort -k1,1 -k4,5n | bgzip > gtf.processed.gz

    vcfanno -p $(getconf _NPROCESSORS_ONLN) -lua /introme/annotations/conf.lua /introme/annotations/gencode.${params.genome_build}.toml $input_vcf | bgzip > ${params.prefix}.variant_info.vcf.gz
    tabix ${params.prefix}.variant_info.vcf.gz

    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(gnomAD_PM_AF<=${params.allele_frequency} || gnomAD_PM_AF='.')" ${params.prefix}.variant_info.vcf.gz | bgzip > ${params.prefix}.variant_info.filtered.vcf.gz
    tabix ${params.prefix}.variant_info.filtered.vcf.gz

    # Prepare files for next step 
    gunzip -k ${params.prefix}.variant_info.filtered.vcf.gz # MMSplice and SpliceAI needs unzipped input files

    # Remove annotations which interfere with score extraction
    bcftools view -h ${params.prefix}.variant_info.filtered.vcf.gz > ${params.prefix}.variant_info.filtered.rmanno.vcf
    grep -v "^#" ${params.prefix}.variant_info.filtered.vcf | awk '{$8="."; print }' OFS='\t' >> ${params.prefix}.variant_info.filtered.rmanno.vcf
    bgzip ${params.prefix}.variant_info.filtered.rmanno.vcf
    """
}