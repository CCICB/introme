process quality_filter {
    container "${params.data_preprocessing_docker_container}"
    beforeScript 'echo Starting quality_filter'
    afterScript  'echo Completed quality_filter'
    publishDir (path: "${params.outdir}/quality_filter")

    input:
        path input_vcf

    output:
        path "${params.prefix}.quality_filter.vcf.gz", emit: quality_filter
        path "${params.prefix}.quality_filter.vcf.gz.tbi", emit: quality_filter_tbi
        val variant_count

    script:
        """
        bcftools filter --threads \$(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL${params.min_QUAL} || QUAL='.') && MAX(FORMAT/DP[*])${params.min_DP} && MAX(FORMAT/AD[*:1])${params.min_AD}" $input_vcf | bgzip > ${params.prefix}.quality_filter.vcf.gz
        tabix -p vcf ${params.prefix}.quality_filter.vcf.gz
        variant_count=\$(bcftools view -H ${params.prefix}.quality_filter.vcf.gz | wc -l | tr -d ' ')
        echo \$(date +%x_%r) 'Quality filtering complete -' \$(variant_count) 'variants remaining'

        if [[ \$(variant_count) == 0 ]]; then
            echo \$(date +%x_%r) 'No variants passed quality filtering - perhaps rerun with -q (triggers no quality filtering)'
            exit 1
        fi
        """
}