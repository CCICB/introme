process pangolin {
    label 'gpu'
    container "${params.pangolin_docker_container}"
    cpus params.pangolin_cpus
    // containerOptions '--workdir /'

    beforeScript 'echo Starting pangolin'
    afterScript  'echo Completed pangolin'
    publishDir (path: "${params.outdir}/pangolin")

    debug  true

    input:
        path vcf
        path ref_genome
    
    output:
        path "${params.prefix}.pangolin.vcf", emit:  pangolin_output
        // path "pangolin.vcf.gz.tbi", emit: pangolin_output_tbi

    script:
        """
        wget https://compbio.ccia.org.au/introme/files/${params.genome_build}/${params.pangolin_db} --no-check-certificate
        pangolin ${vcf} ${ref_genome} ${params.pangolin_db} ${params.prefix}.pangolin
        
        # bgzip pangolin.vcf
        # tabix pangolin.vcf.gz
        """

}