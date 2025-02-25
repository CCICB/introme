process mmsplice {
    label 'gpu'
    container "${params.mmsplice_docker_container}"
    cpus params.mmsplice_cpus
    memory params.mmsplice_mem

    // containerOptions '--workdir /MMSplice_MTSplice'

    beforeScript 'echo Starting mmsplice'
    afterScript  'echo Completed mmsplice'
    publishDir (path: "${params.outdir}/mmsplice")

    debug  true

    input:
        path vcf
        path ref_genome
        path gtf

    output:
        path "${params.prefix}.mmsplice.vcf", emit:  mmsplice_output
        // path "mmsplice.vcf.gz.tbi", emit: mmsplice_output_tbi

    script:
        """
        pwd
        ls

        export TF_FORCE_GPU_ALLOW_GROWTH=true

		run_mmsplice.py --vcf $vcf --fasta $ref_genome --gtf $gtf --output ${params.prefix}.mmsplice.vcf
        # bgzip mmsplice.vcf
        # tabix mmsplice.vcf.gz
        """
}