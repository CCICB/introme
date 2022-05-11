#!/bin/bash

# Hard filter cutoffs for each annotation, filtering thresholds will depend on your application, we urge you to carefully consider these
min_QUAL='>=200' # The QUAL VCF field
min_DP='>=20' # The sample with the highest depth must have a equal or larger depth than this value
min_AD='>=5' # The sample with the highest number of alternate reads must have a equal or larger number than this value

# Define inputs for subsetting portion of introme
while getopts "b:p:q:v:" opt; do
    case $opt in
        b) input_BED="$OPTARG";; # Input BED file (i.e. regions of interest)
        p) prefix="$OPTARG";; # Output file prefix
        q) no_qual_filter="$OPTARG";; # trigger no quality filter command
        v) input_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list


echo $(gzip -d -c $input_VCF | grep -v '^#' | wc -l) 'variants in input VCF'

# Correctly formats the input VCF file for Introme
bcftools sort $input_VCF | uniq | bgzip > $prefix.sorted.vcf.gz # Ensures the file is sorted correctly prior to subsetting
bcftools norm -m-both $prefix.sorted.vcf.gz | bgzip > $prefix.sorted.norm.vcf.gz # Removes multiallelics

echo $(gzip -d -c $prefix.sorted.norm.vcf.gz | grep -v '^#' | wc -l) 'variants prior to subsetting'

echo $(date +%x_%r) 'Beginning subsetting to genomic regions of interest'
bedtools intersect -header -u -a $prefix.sorted.norm.vcf.gz -b $input_BED | bgzip > $prefix.subset.vcf.gz # -u for unique record in VCF

tabix -p vcf $prefix.subset.vcf.gz

echo $(gzip -d -c $prefix.subset.vcf.gz | grep -v '^#' | wc -l) 'variants after subsetting'
echo $(date +%x_%r) 'Subsetting complete'


echo $(date +%x_%r) 'Beginning quality filtering'

if [ $no_qual_filter == 0 ]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.') && MAX(FORMAT/DP[*])$min_DP && MAX(FORMAT/AD[*:1])$min_AD" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
else
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.')" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
fi

tabix -p vcf $prefix.subset.highquality.vcf.gz

echo $(gzip -d -c $prefix.subset.highquality.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Quality filtering complete'
