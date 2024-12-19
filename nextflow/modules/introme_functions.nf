process introme_functions {
    container "${params.introme_functions_docker_container}"
    // containerOptions '--workdir /'
    beforeScript 'echo Starting introme_functions'
    afterScript  'echo Completed introme_functions'
    publishDir (path: "${params.outdir}/introme_functions")

    input:
        path ag_script_path
        path ese_script_path

        path annotated_vcf
        path rmanno_vcf
        path ref_genome
        path template_header_vcf
		 
    output:
        path "${params.prefix}.introme_annotate.ag_check.vcf.gz", emit: ag_check
        path "${params.prefix}.introme_annotate.ag_check.vcf.gz.tbi", emit: ag_check_tbi
        path "${params.prefix}.introme_annotate.ESE.vcf", emit: ese_score

    script:
    """
    echo executing pwd
    pwd
    echo executing ls
    ls

    # ./AG_check.sh -i \$out_dir/working_files/\$prefix.subset.highquality.annotated.vcf.gz -r \$reference_genome

    python3 ${ag_script_path} ${annotated_vcf} ${ref_genome} ${template_header_vcf} ${params.prefix}.introme_annotate.ag_check.vcf
    bgzip ${params.prefix}.introme_annotate.ag_check.vcf
    tabix ${params.prefix}.introme_annotate.ag_check.vcf.gz

    # python ESE_ESS_scoring.py \$out_dir/working_files/\$prefix.subset.highquality.annotated.filtered.rmanno.vcf \$out_dir/working_files/\$prefix.subset.highquality.ESE.vcf \$reference_genome
    ## NOTE here, the previous version of ESE scoring removed the header lines, so next step is to prepend it back
    # cat annotations/introme_annotate.vcf \$out_dir/working_files/\$prefix.subset.highquality.ESE.vcf | bgzip > introme_annotate.ESE.vcf

    python3 ${ese_script_path} ${annotated_vcf} ./${params.prefix}.introme_annotate.ESE.vcf ${ref_genome}
    """

}