process spip {
    container "${params.spip_docker_container}"
    containerOptions "--workdir / -v ${params.outdir}/variant_info/:/data/"
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
        ls
        pwd
        Rscript ./SPiP/SPiPv2.1_main.r -I /data/${vcf} -O spip.vcf -g ${params.genome_build} --VCF
        bgzip spip.vcf
        tabix spip.vcf.gz
        """

}