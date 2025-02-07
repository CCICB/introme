process splicing_anno {
    container "${params.introme_functions_docker_container}"
    beforeScript 'echo Starting splicing_anno'
    afterScript  'echo Completed splicing_anno'
    publishDir (path: "${params.outdir}/splicing_anno")

    debug true

    input:
        path vcf
        path conf_lua
        path ensemble_anno_toml

		// path cadd
		// path cadd_tbi
		// path dbscSNV
		// path dbscSNV_tbi
		// path branchpointer
		// path branchpointer_tbi

        path spliceai_output
		// path spliceai_output_tbi
        path mmsplice_output
		// path mmsplice_output_tbi
        path pangolin_output
        // path pangolin_output_tbi
        path spip_output
        // path spip_output_tbi
        path spliceogen_output
        //path squirl_output
        //path squirl_output_tbi
		path ag_check
		path ag_check_tbi
		path ese_score
		path ese_score_tbi

    output:
		path "${params.prefix}.highquality.annotated.filtered.ensemblescored.vcf.gz", emit: splicing_anno_output
        // path "${params.prefix}.annotated.tsv", emit: annotated_tsv

    script:
    """
    echo splicing_anno
    sed -i -E 's/([0-9]+),(NM_|ENST|ENSG)/\\1\\&\\2/g' ${spliceai_output} #
    bgzip -c ${spliceai_output} > spliceai.vcf.gz
    tabix -p vcf spliceai.vcf.gz

    bgzip -c ${mmsplice_output} > mmsplice.vcf.gz
    tabix -p vcf mmsplice.vcf.gz

    bgzip -c ${pangolin_output} > pangolin.vcf.gz
    tabix -p vcf pangolin.vcf.gz

    # vcfanno needs an = somewhere in header lines
    sed -i -E 's_(##SPiP output) (v[0-9]+.[0-9]*)_\\1=\\2_' ${spip_output}
    bgzip -c ${spip_output} > spip.vcf.gz
    tabix -p vcf spip.vcf.gz

    bgzip -c ${spliceogen_output} > spliceogen.tsv.gz
    tabix -s1 -b2 -e3 spliceogen.tsv.gz

    ln -s ${ag_check} introme_annotate.ag_check.vcf.gz
    ln -s ${ag_check_tbi} introme_annotate.ag_check.vcf.gz.tbi

    ln -s ${ese_score} introme_annotate.ESE.tsv.gz
    ln -s ${ese_score_tbi} introme_annotate.ESE.tsv.gz.tbi
    
    pwd
    ls


    vcfanno \
        -base-path ./ \
        -p \$(getconf _NPROCESSORS_ONLN) \
        -lua ${conf_lua} \
        ${ensemble_anno_toml} \
        ${vcf} > ${params.prefix}.highquality.annotated.filtered.ensemblescored.vcf
    
    bgzip -k ${params.prefix}.highquality.annotated.filtered.ensemblescored.vcf
    """
    // # vcfanno -lua /introme/annotations/conf.lua /introme/annotations/vcfanno_splicing.toml ${vcf} \
    // #    | bgzip > ${params.prefix}.highquality.annotated.filtered.scored_anno.vcf.gz    
    
    // # wget https://github.com/CCICB/introme/blob/master/annotations/U12.${params.genome_build}.bed.gz
    // # wget https://github.com/CCICB/introme/blob/master/annotations/U12.${params.genome_build}.bed.gz.tbi

    // # vcfanno -lua /introme/annotations/conf.lua /introme/annotations/vcfanno_splicing_run.toml ${params.prefix}.highquality.annotated.filtered.scored_anno.vcf.gz | bgzip > ${params.prefix}.highquality.annotated.filtered.scored.vcf.gz
    
    // # java -jar vcftotsv-assembly-0.1.jar --inputFile ${params.prefix}.highquality.annotated.filtered.scored.vcf.gz --outputFile ${params.prefix}.annotated.tsv

    // # Sort by chromosome and coordinate
    // # sort -k1,1n -k2,2n ${params.prefix}.annotated.tsv
}