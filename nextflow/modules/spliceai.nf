process spliceai {
    container "${params.spliceai_docker_container}"
    containerOptions '--workdir /'
    beforeScript 'echo Starting spliceai'
    afterScript  'echo Completed spliceai'
    publishDir (path: "${params.outdir}/spliceai")
    cpus 8 // On my computer (Gab's) the max available is 8 - should be 32

    debug  true

    input:
        path vcf
        path ref_genome
        val distance // 1000: How to add default values? https://github.com/nextflow-io/nextflow/discussions/3714
        val mask // 0
        
    output:
        path "spliceai.vcf.gz", emit:  spliceai_output
        path "spliceai.vcf.gz.tbi", emit: spliceai_output_tbi
    
    script:
        """
        wget https://compbio.ccia.org.au/introme/files/${params.genome_build}/${params.spliceai_db} --no-check-certificate
        ls
        spliceai -I $vcf -O spliceai.vcf -R $ref_genome -A ${params.spliceai_db} -D $distance -M $mask
        # Check if need to sort 
        bgzip spliceai.vcf 
        tabix spliceai.vcf.gz
        """
}