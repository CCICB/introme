process spliceogen {
    container "${params.spliceogen_docker_container}"
    containerOptions = '--entrypoint ""'
    beforeScript 'echo Starting spliceogen'
    afterScript  'echo Completed spliceogen'
    publishDir (path: "${params.outdir}/spliceogen")
    debug  true

    input:
        path vcf
        path ref_genome
        path gtf
    
    output:
        path "${params.prefix}.spliceogen.tsv", emit:  spliceogen_output
    
    script:
        """
        echo pwd
        pwd
        echo ls
        ls

        ORIG_DIR=\$(pwd)

        # cd ../Spliceogen
		# gunzip -c $vcf > ${params.prefix}.spliceogen_input.vcf
		
        cd /Spliceogen
        ./RUN.sh -input \$(realpath "\$ORIG_DIR/${vcf}") \
                 -fasta \$(realpath "\$ORIG_DIR/${ref_genome}") \
                 -gtf \$(realpath "\$ORIG_DIR/${gtf}")

        ls output

        mv output/*_out.txt \$ORIG_DIR/${params.prefix}.spliceogen.tsv
        """
}   