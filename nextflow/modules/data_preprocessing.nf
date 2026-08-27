process data_preprocessing {
    container "${params.data_preprocessing_docker_container}"
    beforeScript 'echo Starting data_preprocessing'
    afterScript  'echo Completed data_preprocessing'
    publishDir (path: "${params.outdir}/data_preprocessing")

    input: 
        path input_vcf
        path ref_genome
        path input_gtf
        path chrRename
    
    output:
        path "${params.prefix}.subset.vcf.gz", emit: preprocessed_output
        path "${params.prefix}.subset.vcf.gz.tbi", emit: preprocessed_output_tbi
        path "sorted.gtf.gz", emit: sorted_gtf

    script:
        """
        CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

        gtf_path=\$(echo input_gtf | sed 's/.gtf.gz$/\.sorted\.gtf\.gz/')

        echo \$(date +%x_%r) \$(bcftools view -H "$input_vcf" | wc -l) 'variants prior to subsetting'

        bcftools annotate --rename-chrs "$chrRename" -Ou "$input_vcf" \
            | bcftools sort -Ou \
            | bcftools view -t "\$CHROMS" -Ou \
            | bcftools norm -f "$ref_genome" -c e -m-both -Ou \
            | bcftools norm -d exact -Ou \
            | bcftools sort -Oz -o ${params.prefix}.sorted.norm.vcf.gz

        if tabix -f -p gff "$input_gtf"; then
            cp -L $input_gtf "\$gtf_path"
        else
            # Need to sort and bgzip the GTF file before tabix indexing
            echo \$(date +%x_%r) 'GTF file not in BGZF format - sorting and bgzipping for tabix indexing'
            # Print header lines (or override grep result with true) then sort the rest.
            (
                gunzip -c "$input_gtf" | grep '^#' || true
                gunzip -c "$input_gtf" | grep -v '^#' |
                    sort -k1,1 -k4,4n -k5,5n -s
            ) | bgzip > "\$gtf_path"

            tabix -f "\$gtf_path"
        fi

        rm "\$gtf_path.tbi"

        if [ -z ${params.bed} ]; then
            echo \$(date +%x_%r) 'No BED file provided - Beginning subsetting to GTF regions'
            bedtools intersect -header -u -a ${params.prefix}.sorted.norm.vcf.gz -b \$gtf_path | bgzip > ${params.prefix}.subset.vcf.gz # -u for unique record in VCF
        else
            echo \$(date +%x_%r) 'BED file provided - Beginning subsetting to genomic regions of interest'
            bedtools intersect -header -u -a ${params.prefix}.sorted.norm.vcf.gz -b ${params.bed} | bgzip > ${params.prefix}.subset.vcf.gz # -u for unique record in VCF
        fi

        tabix -p vcf ${params.prefix}.subset.vcf.gz

        variant_count=\$(bcftools view -H ${params.prefix}.subset.vcf.gz | wc -l | tr -d ' ')
        
        if [ "\$variant_count" -eq "0" ]; then
            if [ -z ${params.bed} ]; then
                echo \$(date +%x_%r) 'No variants were in regions of interest - is your GTF file restricted to certain regions?'
            else
                echo \$(date +%x_%r) 'No variants were in regions of interest - perhaps expand your BED file to more regions of interest'
            fi
            exit 1
        else
            echo \$(date +%x_%r) "\$variant_count variants prior to subsetting"
        fi
        """
}