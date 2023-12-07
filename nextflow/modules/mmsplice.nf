process mmsplice {
    container "${params.mmsplice_docker_container}"
    containerOptions '--workdir /'
    beforeScript 'echo Starting mmsplice'
    afterScript  'echo Completed mmsplice'
    publishDir (path: "${params.outdir}/mmsplice")
    cpus 4
    memory '10 GB'
    debug  true

    input:
        path vcf
        path ref_genome
        path gtf
    
    output:
        path "mmsplice.vcf.gz", emit:  mmsplice_output
        path "mmsplice.vcf.gz.tbi", emit: mmsplice_output_tbi

    script:
        """
        pwd
        ls
        cd ../MMSplice_MTSplice
		python3 run_mmsplice.py --vcf $vcf --fasta $ref_genome --gtf $gtf --output mmsplice.vcf
        bgzip mmsplice.vcf
        tabix mmsplice.vcf.gz
        """

}