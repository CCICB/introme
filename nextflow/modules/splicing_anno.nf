process splicing_anno {
    container "${params.variant_info_docker_container}"
    beforeScript 'echo Starting splicing_anno'
    afterScript  'echo Completed splicing_anno'
    publishDir (path: "${params.outdir}/splicing_anno")

    input:
        path vcf

		// path cadd
		// path cadd_tbi
		// path dbscSNV
		// path dbscSNV_tbi
		// path branchpointer
		// path branchpointer_tbi

        path spliceai_output
		path spliceai_output_tbi
        path mmsplice_output
		path mmsplice_output_tbi
        path pangolin_output
        path pangolin_output_tbi
        path spip_output
        path spip_output_tbi
        //path squirl_output
        //path squirl_output_tbi
		path functions
		path functions_tbi

    output:
		path "${params.prefix}.highquality.annotated.filtered.scored.vcf.gz", emit: splicing_anno_output
        path "${params.prefix}.annotated.tsv", emit: annotated_tsv

    script:
    """
    vcfanno -lua /introme/annotations/conf.lua /introme/annotations/vcfanno_splicing.toml ${vcf} | bgzip > ${params.prefix}.highquality.annotated.filtered.scored_anno.vcf.gz    
    
    wget https://github.com/CCICB/introme/blob/master/annotations/U12.${params.genome_build}.bed.gz
    wget https://github.com/CCICB/introme/blob/master/annotations/U12.${params.genome_build}.bed.gz.tbi

    vcfanno -lua /introme/annotations/conf.lua /introme/annotations/vcfanno_splicing_run.toml ${params.prefix}.highquality.annotated.filtered.scored_anno.vcf.gz | bgzip > ${params.prefix}.highquality.annotated.filtered.scored.vcf.gz
    
    java -jar vcftotsv-assembly-0.1.jar --inputFile ${params.prefix}.highquality.annotated.filtered.scored.vcf.gz --outputFile ${params.prefix}.annotated.tsv

    # Sort by chromosome and coordinate
    sort -k1,1n -k2,2n ${params.prefix}.annotated.tsv
    """
}