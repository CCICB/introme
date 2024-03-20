process introme_functions {
    container "${params.introme_functions_docker_container}"
    containerOptions '--workdir /'
    beforeScript 'echo Starting introme_functions'
    afterScript  'echo Completed introme_functions'
    publishDir (path: "${params.outdir}/introme_functions")

    input:
        path annotated_vcf 
        path ref_genome
		 
    output:
        path "introme_annotate.functions2.vcf.gz", emit: annotate_functions
        path "introme_annotate.functions2.vcf.gz.tbi", emit: annotate_functions_tbi

    script:
    """
    python ../../AG_check ${annotated_vcf} ${ref_genome}
    bgzip introme_annotate.functions2.vcf
    tabix introme_annotate.functions2.vcf.gz
    """

}