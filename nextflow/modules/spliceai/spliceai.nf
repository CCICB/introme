process spliceai {
    label 'gpu'
    container "${params.spliceai_docker_container}"
    cpus params.spliceai_cpus

    beforeScript 'echo Starting spliceai'
    afterScript  'echo Completed spliceai'
    publishDir (path: "${params.outdir}/spliceai")

    debug  true

    input:
        path vcf
        path ref_genome
        val distance // 1000: How to add default values? https://github.com/nextflow-io/nextflow/discussions/3714
        val mask // 0
        
    output:
        // path "spliceai.vcf.gz", emit:  spliceai_output
        // path "spliceai.vcf.gz.tbi", emit: spliceai_output_tbi
        path "${params.prefix}.spliceai.vcf", emit: spliceai_output
    
    script:
        """
        export TF_FORCE_GPU_ALLOW_GROWTH=true
        wget https://compbio.ccia.org.au/introme/files/${params.genome_build}/${params.spliceai_db} --no-check-certificate
        touch ${params.prefix}.spliceai.vcf

        ls -lah
        spliceai -I ${vcf} -O ${params.prefix}.spliceai.vcf -R ${ref_genome} -A ${params.spliceai_db} -D ${distance} -M ${mask} 1>./log

        sed -i -E 's/([0-9]+),(NM_|ENST|ENSG)/\\1\\&\\2/g' ${params.prefix}.spliceai.vcf

        # Check if need to sort 
        # bgzip spliceai.vcf 
        # tabix spliceai.vcf.gz
        """
}