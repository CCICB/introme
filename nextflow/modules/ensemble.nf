process ensemble_infer {
    container "${params.introme_functions_docker_container}"
    beforeScript 'echo Starting Introme ensemble inference'
    afterScript  'echo Completed Introme ensemble inference'
    publishDir (path: "${params.outdir}/ensemble_infer")

    debug true
    cache false

    input:
        path ensemble_score_script_path
        path clf_model_path
        path columns_path
        path splicing_anno_vcf
        path introme_toml
    
    output:
		path "${params.prefix}.introme.predictions.tsv", emit: ensemble_report
        path "${params.prefix}.introme.predictions.vcf.gz", emit: ensemble_vcf

    script:
    """
    python3 ${ensemble_score_script_path} infer \
        --model-path ${clf_model_path} \
        --columns-path ${columns_path} \
        --input-vcf ${splicing_anno_vcf} \
        --output-tsv ${params.prefix}.introme.predictions.tsv

    bgzip -c ${params.prefix}.introme.predictions.tsv > introme.predictions.tsv.gz
    tabix -S1 -s1 -b2 -e2 introme.predictions.tsv.gz

    introme_col=\$(
        zcat introme.predictions.tsv.gz \
        | head -1 \
        | tr '\t' '\n' \
        | awk '\$0 == "introme_score" {print NR; exit}'
    )

    sed "s/__INTROME_COL__/\${introme_col}/" introme.toml > introme.generated.toml

    vcfanno \
        -base-path ./ \
        -p \$(getconf _NPROCESSORS_ONLN) \
        introme.generated.toml \
        ${splicing_anno_vcf} > ${params.prefix}.introme.predictions.vcf

    bgzip ${params.prefix}.introme.predictions.vcf
    """
}

process ensemble_train {
    container "${params.introme_functions_docker_container}"
    beforeScript 'echo Starting Introme ensemble training'
    afterScript  'echo Completed Introme ensemble training'

    publishDir (path: "${params.ml_save_dir}", mode: 'copy', pattern: 'models/*')
    publishDir (path: "${params.ml_log_dir}", mode: 'copy', pattern: 'logs/*')

    debug true
    cache false

    input:
        path ensemble_score_script_path
        path splicing_anno_vcf
        val test_chroms_arg

    output:
        path "models/*", emit: model_artifacts
        path "logs/*", emit: training_logs

    script:
    """
    mkdir -p models logs
    python3 ${ensemble_score_script_path} train \
        --save-dir models \
        --input-vcf ${splicing_anno_vcf} \
        --log-dir logs \
        --run-name ${params.ml_run_name} \
        --test-chroms ${test_chroms_arg}
    """
}
