process introme_functions {
    container "${params.introme_functions_docker_container}"
    // containerOptions '--workdir /'
    beforeScript 'echo Starting introme_functions'
    afterScript  'echo Completed introme_functions'
    publishDir (path: "${params.outdir}/introme_functions")

    input:
        path ag_script_path
        path annotated_vcf 
        path ref_genome
        path template_header_vcf
		 
    output:
        path "introme_annotate.functions2.vcf.gz", emit: annotate_functions
        path "introme_annotate.functions2.vcf.gz.tbi", emit: annotate_functions_tbi

    script:
    """
    echo executing pwd
    pwd
    echo executing ls
    ls

    python ${ag_script_path} ${annotated_vcf} ${ref_genome} ${template_header_vcf}
    bgzip introme_annotate.functions2.vcf
    tabix introme_annotate.functions2.vcf.gz
    """

}