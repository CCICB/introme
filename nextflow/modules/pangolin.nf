process pangolin {
    container "${params.pangolin_docker_container}"
    beforeScript 'echo Starting pangolin'
    afterScript  'echo Completed pangolin'
    publishDir (path: "${params.outdir}/pangolin")
    debug  true

    input:
        path vcf
        path ref_genome
    
    output:
        path "pangolin.vcf.gz", emit:  pangolin_output
        path "pangolin.vcf.gz.tbi", emit: pangolin_output_tbi

    script:
        """
        wget https://compbio.ccia.org.au/introme/files/${params.genome_build}/${params.spliceai_db} --no-check-certificate
        pangolin $vcf $ref_genome ${params.pangolin_db} pangolin.vcf
        bgzip pangolin.vcf
        tabix pangolin.vcf.gz
        """

}