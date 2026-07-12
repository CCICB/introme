/*
    Run via the below steps:
    ./update_params.sh [options] -b <subset.bed.gz> -p <prefix> -r <reference_genome.fa> -v <variants.vcf.gz> -g <annotations.gtf>
    nextflow run https://github.com/CCICB/introme/tree/nextflow/nextflow -params-file params.json -process.echo true
    
    nextflow run https://github.com/CCICB/introme/tree/nextflow/nextflow --vcf './test.vcf' 

    FOR TESTING [as of 06/12]:
    nextflow run /Users/gyounes/Desktop/introme/nextflow/ --vcf test.vcf.gz --ref_genome hg38.fa --gtf output/data_preprocessing/sorted.gtf.gz --prefix 6_12_test --quality_filter false -params-file /Users/gyounes/Desktop/introme/nextflow/params.json

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////                         Progress Log                       ////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    [20/09/2023] Only runs spliceai for now
    [04/10/2023] Contains code to run data_preprocessing+quality_filtering(subsetting) + programs
                            Still need to add in vcf anno process and fix squirls. Need to determine what degree of
                            parallel running for the programs will be used, and how input variables are given with
                            respect to the shell script. Need to go through the docker of each program and double check
                            that it is up to date.
    [25/10/2023] Steps 1,2,3,4 complete. 
                            TODO's: 
                                    - In step 4, need to load the lua/toml files in through the docker
    [01/11/2023] Modified the lua and toml files that would run during splicing_anno - have not yet been tested.
                            Files live in my Downloads folder - not in the Docker yet
                            Need to fix/look at the spip and squirls database creation
    [14/11/2023] splicing_anno process is at a tentative stage. Changes to be made:
                                    - Upload the toml and lua files to git
                                    - Add the toml and lua files into the docker 
                                    - Determine if CADD/dbSnV is going to be uncommented from vcfanno_splicing_run.toml
    [15/11/2023] Testing steps 1 - 6
                            Data input: it seems that we can't pass files in from CLI so instead they would have to be 
                                                    specified in params.
                                                    we could also have default file paths in a github with this and have them set as
                                                    the params value or use fromPath in the main
                            Removed slurm {} from nextflow.config.
                            TODO:
                                    - Consider what needs to go inside the nextflow.config file
    [06/12/2023] Currently in testing:
                                - data_preprocessing works
                                - variant_info works
                                - currently trying to get spip/spliceai/mmsplice/pangolin/squirl working
                                    - current issue: when spip runs the vcf cannot be found: path[1]="/data/6_12_test.variant_info.filtered.rmanno.vcf": No such file or directory
                                                                     This is because the process is running INSIDE the docker container. This means that there we need to mount the directory
                                                                     containing the output files of the previous nextflow step into the container, so that SPIP can run on the relevant vcf file.
                                                                     The process also has to run inside the docker container, otherwise the spip repo can't be found.
                                                                     Attempted fix:
                                                                        - Adding "-v ${params.outdir}/variant_info/:/data/" to the containerOptions. The intenion of this is to insert the output vcf
                                                                            file into the spip docker container.
                                                                        - Currently throws this error: path[1]="/data/6_12_test.variant_info.filtered.rmanno.vcf": No such file or directory
                                - commented out splicing_anno and introme_functions to focus on debugging the parallel programs 
*/


/*
 * Provide workflow description and default param values to user
 */
def intromeBanner(workflowParams) {
    return """\

====================================================================================
██╗███╗   ██╗████████╗██████╗  ██████╗ ███╗   ███╗███████╗    ██████╗     ██████╗ 
██║████╗  ██║╚══██╔══╝██╔══██╗██╔═══██╗████╗ ████║██╔════╝    ╚════██╗   ██╔═████╗
██║██╔██╗ ██║   ██║   ██████╔╝██║   ██║██╔████╔██║█████╗       █████╔╝   ██║██╔██║
██║██║╚██╗██║   ██║   ██╔══██╗██║   ██║██║╚██╔╝██║██╔══╝      ██╔═══╝    ████╔╝██║
██║██║ ╚████║   ██║   ██║  ██║╚██████╔╝██║ ╚═╝ ██║███████╗    ███████╗██╗╚██████╔╝
╚═╝╚═╝  ╚═══╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝    ╚══════╝╚═╝ ╚═════╝                        
====================================================================================


Runs the Introme pipeline in the following steps:


    1. Data Pre-Processing: Subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)
    2. Quality Filtering: Hard filtering on variant quality (this is here to reduce the number of variants going into the CPU-costly annotation step below)
    3. Annotate the VCF: Annotate the subsetted VCF with useful information, to be used for filtering downstream
    4. Frequency Filtering: Hard filtering on the values of annotations added in the previous step
    5. Execute Programs: Run MMSplice, Splice AI, Pangolin, Spliceogen, Squirls and Spip
    6. Splicing Annotations: Runs annotations on splicing events on the output of all the programs from step 5
    7. Machine Learning: Generate consensus scores using a ML algorithm

Inputs:
    vcf               : ${workflowParams.vcf}
    reference genome  : ${workflowParams.ref_genome}
    gtf               : ${workflowParams.gtf}
    genome build      : ${workflowParams.genome_build}
    prefix            : ${workflowParams.prefix}
    bed               : ${workflowParams.bed}

Process:
        Dockers:
                spliceai            : ${workflowParams.spliceai_docker_container}
                mmsplice            : ${workflowParams.mmsplice_docker_container}
                spliceogen          : ${workflowParams.spliceogen_docker_container}
                pangolin            : ${workflowParams.pangolin_docker_container}
                spip                : ${workflowParams.spip_docker_container}
                squirl              : ${workflowParams.squirl_docker_container}
                data_preprocessing  : ${workflowParams.data_preprocessing_docker_container}
                variant_info        : ${workflowParams.variant_info_docker_container}
                introme_functions   : ${workflowParams.introme_functions_docker_container}

Output:
        Output folder  : ${workflowParams.outdir}
        
"""
}

def resolveRequiredPath(workflowParams, key) {
    def value = workflowParams[key]
    if (value == null || value.toString().trim() == '') {
        error "Missing required parameter --${key}"
    }
    file(value.toString(), checkIfExists: true)
}

def resolveRequiredPatternPath(workflowParams, key, token) {
    def pattern = workflowParams[key]
    if (pattern == null || pattern.toString().trim() == '') {
        error "Missing required parameter --${key}"
    }
    def resolved = String.format(pattern.toString(), token)
    file(resolved, checkIfExists: true)
}


/* 
 * Import modules 
 */
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { data_preprocessing }  from './modules/data_preprocessing.nf'
include { quality_filter }      from './modules/quality_filter.nf'
include { variant_info }        from './modules/variant_info.nf'
include { spliceai }            from './modules/spliceai/spliceai.nf'
include { mmsplice }            from './modules/mmsplice/mmsplice.nf'
include { spliceogen }          from './modules/spliceogen/spliceogen.nf'
include { pangolin }            from './modules/pangolin/pangolin.nf'
include { spip }                from './modules/spip/spip.nf'
include { squirl }              from './modules/squirl.nf'
include { introme_functions }   from './modules/introme_functions.nf'
include { splicing_anno }       from './modules/splicing_anno.nf'
include { ensemble_infer; ensemble_train } from './modules/ensemble.nf'

/* 
 * Main pipeline logic
 */
workflow {
    log.info intromeBanner(params)
    log.info paramsSummaryLog(workflow)

    // Validate required top-level params once at workflow start.
    def ml_mode = (params.ml_mode ?: 'infer').toString().toLowerCase()
    def ml_test_chroms_arg = null

    if (!(ml_mode in ['infer', 'train'])) {
        error "Invalid --ml_mode '${params.ml_mode}'. Supported values: infer | train"
    }

    if (ml_mode == 'train') {
        def ml_test_chroms = params.ml_test_chroms
        if (!(ml_test_chroms instanceof List)) {
            ml_test_chroms = ml_test_chroms?.toString()?.split(/[\s,]+/)?.findAll { chrom -> chrom } ?: []
        }

        ml_test_chroms = ml_test_chroms.collect { chrom -> chrom.toString().trim() }.findAll { chrom -> chrom }
        if (ml_test_chroms.isEmpty()) {
            error "When --ml_mode train, --ml_test_chroms must contain at least one chromosome"
        }

        def validChromPattern = ~/^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$/
        def invalidChroms = ml_test_chroms.findAll { chrom -> !(chrom ==~ validChromPattern) }
        if (!invalidChroms.isEmpty()) {
            error "Invalid --ml_test_chroms values: ${invalidChroms.join(', ')}. Allowed values are chr1-chr22, chrX, chrY"
        }

        ml_test_chroms_arg = ml_test_chroms.join(' ')
    }

    // Input variables
    def vcf = resolveRequiredPath(params, 'vcf')
    def ref_genome = resolveRequiredPath(params, 'ref_genome')
    def gtf = resolveRequiredPath(params, 'gtf')
    def chrRename = resolveRequiredPath(params, 'chrRename')

    // STEP 1: subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)
    data_preprocessing(vcf, ref_genome, gtf, chrRename)

    // STEP 2: Hard filtering on variant quality (this is here to reduce the number of variants going into the CPU-costly annotation step below)
    def anno_input
    def run_quality_filter = params.quality_filter.toString().toBoolean()
    if (run_quality_filter) {
        // run filter process
        quality_filter(data_preprocessing.out.preprocessed_output)
        // trigger input file for next step to be filtered output
        anno_input = quality_filter.out.quality_filter
    } else {
        // Note: might need to move this outside of the if statement
        // input file for next step is output of data_preprocessing
        anno_input = data_preprocessing.out.preprocessed_output
    }

    // STEP 3: annotate the subsetted VCF with useful information, to be used for filtering downstream
    //         and run hard filtering on the values of annotations added in the previous step
    def conf_pre_lua_path = resolveRequiredPath(params, 'conf_pre_anno_lua')
    def toml_path = resolveRequiredPatternPath(params, 'gencode_toml_pattern', params.genome_build.toString())
    variant_info(anno_input, data_preprocessing.out.sorted_gtf, conf_pre_lua_path, toml_path)


    // STEP 4: Run MMSplice, Splice AI, Pangolin, Spliceogen, Squirls and Spip

    // Define paramaters for SpliceAI
    def distance = 1000
    def mask = 0
    // Run SpliceAI
    spliceai(variant_info.out.variant_info_rmanno, ref_genome, distance, mask)

    // // Run MMSplice
    mmsplice(variant_info.out.variant_info_rmanno, ref_genome, gtf)

    // Run Pangolin
    pangolin(variant_info.out.variant_info_rmanno, ref_genome)

    // Run Spip
    spip(variant_info.out.variant_info_rmanno)

    // Run Squirl
    // download from patricia server to run squirl??? 
    // TODO fix squirl
    // SQUIRLS_DATA = resolveRequiredPath(params, 'squirls_data_dir')
    // squirl(SQUIRLS_DATA, variant_info.out.variant_info_rmanno)

    // Run Splicoegen
    spliceogen(variant_info.out.variant_info_rmanno, ref_genome, gtf)

    // STEP 5: Execute introme functions such as AG_check
    def ag_script_path = resolveRequiredPath(params, 'ag_script_path')
    def ese_script_path = resolveRequiredPath(params, 'ese_script_path')
    // mnv_script_path = file('../MNV.sh')
    def template_header_vcf = resolveRequiredPath(params, 'template_header_vcf')

    introme_functions(ag_script_path, ese_script_path,
                                        // variant_info.out.variant_info,
                                        variant_info.out.variant_info_stripped,
                                        // variant_info.out.variant_info_rmanno,
                                        ref_genome,
                                        template_header_vcf)

    def conf_ensemble_lua_path = resolveRequiredPath(params, 'conf_ensemble_lua')
    def ensemble_anno_toml = resolveRequiredPath(params, 'ensemble_anno_toml')
    def annotate_toml = resolveRequiredPatternPath(params, 'annotate_toml_pattern', params.genome_build.toString())

    def branchpointer_dir = resolveRequiredPath(params, 'branchpointer_dir')
    def regions_dir = resolveRequiredPath(params, 'regions_dir')
    def u12_dir = resolveRequiredPath(params, 'u12_dir')
    
    // STEP 6: Run splicing annotations
    splicing_anno(
        variant_info.out.variant_info, // .vcf.gz
        conf_ensemble_lua_path,
        annotate_toml,
        ensemble_anno_toml,

        branchpointer_dir,
        regions_dir,
        u12_dir,

        // tbi files generated within splicing_anno module
        spliceai.out.spliceai_output, 
        mmsplice.out.mmsplice_output,
        pangolin.out.pangolin_output,
        spip.out.spip_output,
        spliceogen.out.spliceogen_output,
        //squirl.out.squirl_output,

        introme_functions.out.ag_check,
        introme_functions.out.ag_check_tbi,
        introme_functions.out.ese_score,
        introme_functions.out.ese_score_tbi
    )

    // STEP 7: Generate consensus scores - ML mode selection (infer | train)
    def ensemble_score_script_path = resolveRequiredPath(params, 'ensemble_score_script_path')

    if (ml_mode == 'infer') {
        def clf_model_path = resolveRequiredPath(params, 'ml_model_path')
        def columns_path = resolveRequiredPath(params, 'ml_columns_path')

        ensemble_infer(
            ensemble_score_script_path,
            clf_model_path,
            columns_path,
            splicing_anno.out.splicing_anno_output
        )
    } else {
        ensemble_train(
            ensemble_score_script_path,
            splicing_anno.out.splicing_anno_output,
            ml_test_chroms_arg
        )
    }
}
