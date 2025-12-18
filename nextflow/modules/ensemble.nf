process ensemble {
    container "${params.introme_functions_docker_container}"
    beforeScript 'echo Starting Introme ensemble score generation'
    afterScript  'echo Completed Introme ensemble score generation'
    publishDir (path: "${params.outdir}/ensemble")

    debug true

    input:
        path ensemble_score_script_path
        path clf_model_path
        path columns_path
        path splicing_anno_vcf
    
    output:
		path "${params.prefix}.introme.predictions.tsv", emit: ensemble_output

    script:
    """
    python3 ${ensemble_score_script_path} ${clf_model_path} ${columns_path} ${splicing_anno_vcf} ${params.prefix}.introme.predictions.tsv
    """
}
