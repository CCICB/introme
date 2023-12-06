process spip {
    container "${params.spip_docker_container}"
    beforeScript 'echo Starting spip'
    afterScript  'echo Completed spip'
    publishDir (path: "${params.outdir}/spip")
    debug  true

    input:
        path vcf

    output:
        path "spip.vcf.gz", emit:  spip_output
        path "spip.vcf.gz.tbi", emit: spip_output_tbi

    script:
        """
        pwd
        ls
        Rscript ./SPiP/SPiPv2.1_main.r -I $vcf -O spip.vcf -g ${params.genome_build} --VCF
        bgzip spip.vcf
        tabix spip.vcf.gz
        """

}