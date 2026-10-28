process pangolin {
    label 'gpu'
    container "${params.pangolin_docker_container}"
    cpus params.pangolin_cpus
    memory params.pangolin_mem
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

    script:
        """
        # TODO: add generically named softlink to db and refer to that so that the gencode version is not hardcoded in the command
        pangolin-batched --batch_size 8 -d 1000 ${vcf} ${ref_genome} /opt/introme/db/gencode.v50.annotation.all-transcripts.pangolin.db ${params.prefix}.pangolin
        
        # bgzip pangolin.vcf
        # tabix pangolin.vcf.gz
        """
}