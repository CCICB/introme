/*
  Run via the below steps:
  ./update_params.sh [options] -b <subset.bed.gz> -p <prefix> -r <reference_genome.fa> -v <variants.vcf.gz>"
  nextflow run https://github.com/CCICB/introme/tree/nextflow/nextflow -params-file params.json -process.echo true

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
*/



/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Provide workflow description and default param values to user
 */
log.info """\

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
  vcf               : ${params.vcf}
  reference genome  : ${params.ref_genome} 
  gtf               : ${params.gtf}
  genome build      : ${params.genome_build}
  prefix            : ${params.prefix}

Process:
    Dockers:
        spliceai            : ${params.spliceai_docker_container}
        mmsplice            : ${params.mmsplice_docker_container}
        spliceogen          : ${params.spliceogen_docker_container}
        pangolin            : ${params.pangolin_docker_container}
        spip                : ${params.spip_docker_container}
        squirl              : ${params.squirl_docker_container}
        data_preprocessing  : ${params.data_preprocessing_docker_container}
        variant_info        : ${params.variant_info_docker_container}
        introme_functions   : ${params.introme_functions_docker_container}

Output:
    Output folder  : ${params.outdir}
    
"""


/* 
 * Import modules 
 */
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { data_preprocessing } from './modules/data_preprocessing.nf'
include { quality_filter } from './modules/quality_filter.nf'
include { variant_info } from './modules/variant_info.nf'
include { spliceai } from './modules/spliceai.nf'
include { mmsplice } from './modules/mmsplice.nf'
include { spliceogen } from './modules/spliceogen.nf'
include { pangolin } from './modules/pangolin.nf'
include { spip } from './modules/spip.nf'
include { squirl } from './modules/squirl.nf'
include { introme_functions } from './modules/introme_functions.nf'
include { splicing_anno } from './modules/splicing_anno.nf'

/* 
 * Print summary of supplied parameters
 */
log.info paramsSummaryLog(workflow)


/* 
 * Main pipeline logic
 */
workflow {

    // Input variables
    vcf = Channel.fromPath(params.vcf) 
    ref_genome = Channel.fromPath(params.ref_genome) 
    gtf = Channel.fromPath(params.gtf) 
    bed = Channel.fromPath(params.bed)

    // TODO: workout if we need to grab the below
    // File gnomad
    // File gnomad_tbi
    // File regions
    // File regions_tbi
    // File cadd
    // File cadd_tbi
    // File dbscSNV
    // File dbscSNV_tbi
    // File branchpointer
    // File branchpointer_tbi
    // File gencode
    // File gencode_tbi

    // STEP 1: subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)
    data_preprocessing(vcf.first(), bed.first(), ref_genome.first(), gtf.first())

    // STEP 2: Hard filtering on variant quality (this is here to reduce the number of variants going into the CPU-costly annotation step below)
    if (params.quality_filter == true) {
      // run filter process
      quality_filter(vcf.first())
      // trigger input file for next step to be filtered output
      anno_input = quality_filter.out.quality_filter
    } else {
      // Note: might need to move this outside of the if statement
      // input file for next step is output of data_preprocessing
      anno_input = data_preprocessing.out.preprocessed_output
    }

    // STEP 3: annotate the subsetted VCF with useful information, to be used for filtering downstream
    //         and run hard filtering on the values of annotations added in the previous step
    variant_info(anno_input, gtf.first())


    // STEP 4: Run MMSplice, Splice AI, Pangolin, Spliceogen, Squirls and Spip

    // Define paramaters for SpliceAI
    distance = 1000
    mask = 0
    // Run SplicAI
    spliceai(variant_info.out.variant_info_rmanno, ref_genome.first(), distance, mask)

    // Run MMSplice
    mmsplice(variant_info.out.variant_info_rmanno, ref_genome.first(), gtf.first())

    // Run Pangolin
    pangolin(variant_info.out.variant_info_rmanno, ref_genome.first())

    // Run Spip
    spip(variant_info.out.variant_info_rmanno)

    // Run Squirl
    // download from patricia server to run squirl??? 
    // TODO fix squirl
    squirl(SQUIRLS_DATA, variant_info.out.variant_info_rmanno)

    // Run Splicoegen
    spliceogen(variant_info.out.variant_info_rmanno, ref_genome.first(), gtf.first())

    // STEP 5: Execute introme functions such as AG_check
    introme_functions(variant_info.out.variant_info, ref_genome.first())

    // STEP 6: Run splicing annotations
    splicing_anno(
      variant_info.out.variant_info, 
      spliceai.out.spliceai_output, 
      spliceai.out.spliceai_output_tbi,
      mmsplice.out.mmsplice_output,
      mmsplice.out.mmsplice_output_tbi,
      pangolin.out.pangolin_output,
      pangolin.out.pangolin_output_tbi,
      spip.out.spip_output,
      spip.out.spip_output_tbi,
      squirl.out.squirl_output_tbi,
      introme_functions.out.annotate_functions,
      introme_functions.out.annotate_functions_tbi
    )

    // STEP 7: Generate consensus scores - ML
}
