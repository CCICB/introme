#!/bin/bash

# Define inputs for subsetting portion of introme
while getopts "b:p:v:" opt; do
    case $opt in
        b) input_BED="$OPTARG";; # Input BED file (i.e. regions of interest)
        p) prefix="$OPTARG";; # Output file prefix
        v) input_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list


echo $(gzip -d -c $input_VCF | grep -v '^#' | wc -l) 'variants in input VCF'

# Correctly formats the input VCF file for Introme
bcftools sort $input_VCF | uniq | bgzip > $prefix.sorted.vcf.gz # Ensures the file is sorted correctly prior to subsetting
bcftools norm -m-both $prefix.sorted.vcf.gz | bgzip > $prefix.sorted.norm.vcf.gz # Removes multiallelics

echo $(gzip -d -c $input_VCF | grep -v '^#' | wc -l) 'variants prior to subsetting'

echo $(date +%x_%r) 'Beginning subsetting to genomic regions of interest'
bedtools intersect -header -u -a $prefix.sorted.norm.vcf.gz -b $input_BED | bgzip > $prefix.subset.vcf.gz # -u for unique record in VCF

tabix -p vcf $prefix.subset.vcf.gz

echo $(gzip -d -c $prefix.subset.vcf.gz | grep -v '^#' | wc -l) 'variants after subsetting'
echo $(date +%x_%r) 'Subsetting complete'
