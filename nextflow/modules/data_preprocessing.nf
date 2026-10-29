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
        val gtf_prefix
    
    output:
        path "${params.prefix}.subset.vcf.gz", emit: preprocessed_output
        path "${params.prefix}.subset.vcf.gz.tbi", emit: preprocessed_output_tbi
        path "${gtf_prefix}.sorted.gtf.gz", emit: sorted_gtf
        path "${gtf_prefix}.subsetted.gtf.gz", emit: subsetted_gtf

    script:
        def sorted_gtf_name = "${gtf_prefix}.sorted.gtf.gz"

        """
        set -o pipefail
        CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

        gtf_path="${sorted_gtf_name}"

        echo \$(date +%x_%r) \$(bcftools view -H "$input_vcf" | wc -l) 'variants prior to subsetting'

        # Sort and norm the VCF file, and subset to the chromosomes of interest.
        bcftools annotate --rename-chrs "$chrRename" -Ou "$input_vcf" \
            | bcftools sort -Ou \
            | bcftools view -t "\$CHROMS" -Ou \
            | bcftools norm -f "$ref_genome" -c e -m-both -Ou \
            | bcftools norm -d exact -Ou \
            | bcftools sort -Oz -o ${params.prefix}.sorted.norm.vcf.gz

        echo "\$(date +%x_%r) Starting GTF preparation"
        # Sort and index the GTF file.
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

        tabix -f -p gff "\$gtf_path"
        echo "\$(date +%T) Finished GTF preparation"

        # Strategy for small VCF inputs
        # Extract all genomic intervals occupied by variants
        # and retrieve only GTF gene records overlapping those intervals.
        echo \$(date +%x_%r) 'Retrieving GTF records overlapping with variants'
        bcftools query -f '%CHROM\t%POS\t%END\n' \
            "${params.prefix}.sorted.norm.vcf.gz" \
            > variant_regions.tsv
        tabix -R variant_regions.tsv "\$gtf_path" \
            | awk -F '\t' '\$3 == "gene"' \
            | sort -u \
            > hit_genes.gtf
        echo "Input variants are across \$(wc -l < hit_genes.gtf) genes"

        if [ -z ${params.bed} ]; then
            echo \$(date +%x_%r) 'No BED file provided - Beginning subsetting to GTF regions'
            bedtools intersect -header -u -a ${params.prefix}.sorted.norm.vcf.gz -b hit_genes.gtf | bgzip > ${params.prefix}.subset.vcf.gz # -u for unique record in VCF
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

        subset_gtf() {
            local VCF="\$1"
            local GTF="\$2"
            local OUT="\$3"

            # 1. Extract variant intervals from the FINAL subsetted VCF.
            bcftools query -f '%CHROM\\t%POS\\t%END\\n' "\$VCF" \
                > selected_variant_regions.tsv

            # 2. Retrieve gene records overlapping those variants.
            tabix -R selected_variant_regions.tsv "\$GTF" \
                | awk -F '\\t' '\$3 == "gene"' \
                | sort -u \
                > selected_genes.gtf

            # 3. Extract unique gene IDs and whole-gene intervals.
            : > selected_gene_ids.txt
            : > selected_gene_regions.tsv

            awk -F '\\t' '
                match(\$9, /gene_id "[^"]+"/) {
                    id = substr(\$9, RSTART + 9, RLENGTH - 10)

                    if (!seen[id]++) {
                        print id > "selected_gene_ids.txt"
                        print \$1 "\\t" \$4 "\\t" \$5 \
                            > "selected_gene_regions.tsv"
                    }
                }
            ' selected_genes.gtf

            echo "\$(date +%x_%r) Selected \$(wc -l < selected_gene_ids.txt) genes"

            # 4. Retrieve all features in the selected gene regions.
            # Keep only records belonging to selected gene IDs.
            {
                tabix -H "\$GTF"

                if [ -s selected_gene_regions.tsv ]; then
                    tabix -R selected_gene_regions.tsv "\$GTF" \
                        | awk -F '\\t' '
                            NR == FNR {
                                keep[\$1] = 1
                                next
                            }

                            match(\$9, /gene_id "[^"]+"/) {
                                id = substr(\$9, RSTART + 9, RLENGTH - 10)

                                if (id in keep && !seen[\$0]++)
                                    print
                            }
                        ' selected_gene_ids.txt - \
                        | LC_ALL=C sort -s -k1,1 -k4,4n -k5,5n
                fi
            } | bgzip -c > "\$OUT"
        }

        time subset_gtf "${params.prefix}.subset.vcf.gz" "\$gtf_path" "${gtf_prefix}.subsetted.gtf.gz"
        """
}