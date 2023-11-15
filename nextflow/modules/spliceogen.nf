process spliceogen {
    container "${params.spliceogen_docker_container}"
    beforeScript 'echo Starting spliceogen'
    afterScript  'echo Completed spliceogen'
    publishDir (path: "${params.outdir}/spliceogen")
    debug  true

    input:
        path vcf
        path ref_genome
        path gtf
    
    output:
        path "${params.prefix}.spliceogen.txt", emit:  spliceogen_output
    
    script:
        """
        cd ../Spliceogen
		gunzip -c $vcf > ${params.prefix}.spliceogen_input.vcf
		./RUN.sh -input ${params.prefix}.spliceogen_input.vcf -fasta $ref_genome -gtf $gtf
        """

}   