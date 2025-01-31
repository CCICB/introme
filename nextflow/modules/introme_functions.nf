process introme_functions {
    container "${params.introme_functions_docker_container}"
    // containerOptions '--workdir /'
    beforeScript 'echo Starting introme_functions'
    afterScript  'echo Completed introme_functions'
    publishDir (path: "${params.outdir}/introme_functions")

    debug true
    input:
        path ag_script_path
        path ese_script_path

        path variant_info
        path variant_info_stripped
        path variant_info_rmanno
        path ref_genome
        path template_header_vcf
		 
    output:
        path "${params.prefix}.introme_annotate.ag_check.vcf.gz", emit: ag_check
        path "${params.prefix}.introme_annotate.ag_check.vcf.gz.tbi", emit: ag_check_tbi
        path "${params.prefix}.introme_annotate.ESE.tsv.gz", emit: ese_score
        path "${params.prefix}.introme_annotate.ESE.tsv.gz.tbi", emit: ese_score_tbi

    script:
    """
    #####
    # run_introme.sh (step 6)
    echo "Current directory: \$(pwd)"
    echo executing ls
    ls

    ### AG_check ###

    # ./AG_check.sh -i \$out_dir/working_files/\$prefix.subset.highquality.annotated.vcf.gz -r \$reference_genome

    python3 ${ag_script_path} ${variant_info_stripped} ${ref_genome} ${template_header_vcf} ${params.prefix}.introme_annotate.ag_check.vcf
    bgzip -k ${params.prefix}.introme_annotate.ag_check.vcf
    tabix -p vcf ${params.prefix}.introme_annotate.ag_check.vcf.gz

    ### ESE scoring ###

    # python ESE_ESS_scoring.py \$out_dir/working_files/\$prefix.subset.highquality.annotated.filtered.rmanno.vcf \$out_dir/working_files/\$prefix.subset.highquality.ESE.vcf \$reference_genome
    ## NOTE here, the previous version of ESE scoring removed the header lines, so next step is to prepend it back
    # cat annotations/introme_annotate.vcf \$out_dir/working_files/\$prefix.subset.highquality.ESE.vcf | bgzip > introme_annotate.ESE.vcf

    python3 ${ese_script_path} ${variant_info_stripped} ./${params.prefix}.introme_annotate.ESE.vcf ${ref_genome} ${template_header_vcf}
    mv ${params.prefix}.introme_annotate.ESE.vcf ${params.prefix}.introme_annotate.ESE.tsv
    bgzip -k ${params.prefix}.introme_annotate.ESE.tsv
    tabix -s1 -b2 -e2 ${params.prefix}.introme_annotate.ESE.tsv.gz
    """
    // ### MNV scoring ###

    // # Note A:
    // # since there is no subset.highquality.annotated.filtered.vcf.gz (see variant_info.nf),
    // # just the equiv of subset.highquality.annotated.vcf.gz,
    // # ==> which is in nextflow as ${params.prefix}.variant_info.vcf.gz / emitted as variant_info.out.variant_info
    // # MNVs=\$(bcftools filter -i"TYPE!='snp' && TYPE!='indel'" \$out_dir/working_files/\$prefix.subset.highquality.annotated.filtered.vcf.gz | grep -v "^#" | wc -l | tr -d ' ')
    // # MNVs=\$(bcftools filter -i"TYPE!='snp' && TYPE!='indel'" ${variant_info_stripped} | grep -v "^#" | wc -l | tr -d ' ')
    
    // # echo \$MNVs

    // # if [[ \$MNVs > 0 ]]; then
    // #     # Untested as of Dec 18 2024
    // #     echo \$MNVs 'MNV/insdel variants to score'
    // #     echo \${mnv_script_path}
    // #     ls -lah \${mnv_script_path}
    // #     bash \${mnv_script_path} -a ${params.genome_build} -r ${ref_genome} -p ${params.prefix} -f ${variant_info_stripped} #2>/dev/null
    // # else
    // #     echo 'No MNV/insdel variants to score'
    // # fi

    // ### vcfanno ###

    // # vcfanno -p \$(getconf _NPROCESSORS_ONLN) -lua conf.lua annotations/annotate.${params.genome_build}.toml $variant_info 2>/dev/null \
    // #     | bgzip > ${params.prefix}.splicing_anno.vcf.gz
    // # vcfanno -base-path ./ -p \$(getconf _NPROCESSORS_ONLN) -lua ${conf_lua_path} ${annotate_toml_path} $variant_info \
    // #     | bgzip > ${params.prefix}.introme_annotate.splicing_anno.vcf.gz
    // # tabix -f ${params.prefix}.introme_annotate.splicing_anno.vcf.gz

    // ####
    
    // vcfanno -p $(getconf _NPROCESSORS_ONLN) -lua conf.lua annotations/conf_introme.toml $out_dir/working_files/$prefix.subset.highquality.annotated.splicing_anno.vcf.gz 2>/dev/null | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz
    // tabix -f $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz

}