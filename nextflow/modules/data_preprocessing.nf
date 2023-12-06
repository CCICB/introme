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
        echo \$(gzip -d -c $input_vcf | grep -v '^#' | wc -l) 'variants prior to subsetting'

        bcftools annotate --rename-chrs $chrRename $input_vcf | bgzip > ${params.prefix}.chr.vcf.gz
        bcftools sort ${params.prefix}.chr.vcf.gz | bgzip > ${params.prefix}.chr.sort.vcf.gz
        tabix ${params.prefix}.chr.sort.vcf.gz
        bcftools filter -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY ${params.prefix}.chr.sort.vcf.gz | bgzip > ${params.prefix}.chr_filtered.vcf.gz

        bcftools sort ${params.prefix}.chr_filtered.vcf.gz | uniq | bgzip > ${params.prefix}.sorted.vcf.gz  # Ensures the file is sorted correctly prior to subsetting
        bcftools norm -f $ref_genome -c x -m-both ${params.prefix}.sorted.vcf.gz | bgzip > ${params.prefix}.sorted.norm.vcf.gz


        tabix -f $input_gtf && gtf_sorted=1 || gtf_sorted=0

        gtf_path=$input_gtf
        if [ "\$gtf_sorted" -eq "0" ]; then
            gunzip -c $input_gtf | awk '(NR>5)' | sort -k1,1 -k4,4n -k5,5n -s | bgzip > sorted.gtf.gz
            tabix -f sorted.gtf.gz
            gtf_path="sorted.gtf.gz"            
        fi

        echo \$gtf_path

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
        fi
        """
}