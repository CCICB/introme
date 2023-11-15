process squirl {
    container "${params.squirl_docker_container}"
    beforeScript 'echo Starting squirl'
    afterScript  'echo Completed squirl'
    publishDir (path: "${params.outdir}/squirl")
    debug  true

    input:
        path SQUIRLS_DATA // This needs to be downloaded - update dockerfile to accomadate for this
        path vcf

    output:
        path "squirl.vcf.gz", emit:  squirl_output
        path "squirl.vcf.gz.tbi", emit: squirl_output_tbi

    script:
        """
        java -jar squirls-cli/target/squirls-cli-2.0.1.jar annotate-vcf -d $SQUIRLS_DATA $vcf squirl.vcf -f vcf  
        bgzip squirl.vcf
        tabix squirl.vcf.gz
        """
}