process squirl {
    container "${params.squirl_docker_container}"
    // containerOptions '--workdir /Squirls'
    beforeScript 'echo Starting squirl'
    afterScript  'echo Completed squirl'
    publishDir (path: "${params.outdir}/squirl")
    debug  true

    input:
        path SQUIRLS_DATA // This needs to be downloaded on the computer that runs main.nf
        path vcf

    output:
        path "squirl.tsv", emit:  squirl_tsv_output

    script:
        """
        pwd
        ls
        echo $SQUIRLS_DATA
        ORIG_DIR=\$(pwd)

        # cd /Squirls

        # java -jar squirls-cli/target/squirls-cli-2.0.1.jar --help
        java -jar /Squirls/squirls-cli/target/squirls-cli-2.0.1.jar annotate-vcf \
            --report-features -d ${SQUIRLS_DATA} -f tsv ${vcf} squirl 
        """
}