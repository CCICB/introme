#!/bin/bash

# Define inputs for filtering portion of introme
while getopts "a:p:v:" opt; do
    case $opt in
        a) allele_freq="$OPTARG";; # Allele frequency for filtering
        p) prefix="$OPTARG";; # Output file prefix
        v) input_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

echo $(date +%x_%r) 'Beginning filtering by annotation values'

if [ $pop_filter < 1 ]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(MGRB_AF<=$allele_freq || MGRB_AF='.') && (gnomAD_PM_AF<=$allele_freq || gnomAD_PM_AF='.')" $prefix.subset.highquality.annotated.vcf.gz | bgzip > $prefix.subset.highquality.annotated.filtered.vcf.gz
else
    cp $prefix.subset.highquality.annotated.vcf.gz $prefix.subset.highquality.annotated.filtered.vcf.gz
fi
tabix -p vcf $prefix.subset.highquality.annotated.filtered.vcf.gz


echo $(gzip -d -c $prefix.subset.highquality.annotated.filtered.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Frequency filtering complete'
