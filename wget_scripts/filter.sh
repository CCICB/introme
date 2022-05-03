#!/bin/bash

no_pop_filter=0
no_qual_filter=0

# Hard filter cutoffs for each annotation, filtering thresholds will depend on your application, we urge you to carefully consider these
MGRB_AF='<=0.01' # Minimum MGRB allele frequency (healthy Australian population)
gnomad_popmax_AF='<=0.01' # Minimum gnomAD allele frequency in the population with the highest allele frequency in gnomAD
CADD_Phred='>=0' # Minimum CADD score (phred-scaled)
min_QUAL='>=200' # The QUAL VCF field
min_DP='>=20' # The sample with the highest depth must have a equal or larger depth than this value
min_AD='>=5' # The sample with the highest number of alternate reads must have a equal or larger number than this value


# Define inputs for filtering portion of introme
while getopts "fp:qv:" opt; do
    case $opt in
        f) no_pop_filter=1;; # trigger no population filter command
        p) prefix="$OPTARG";; # Output file prefix
        q) no_qual_filter=1;; # trigger no quality filter command
        v) input_VCF="$OPTARG";; # Input VCF file
        
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

echo $(date +%x_%r) 'Beginning quality filtering'

if [ $no_qual_filter == 0 ]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.') && MAX(FORMAT/DP[*])$min_DP && MAX(FORMAT/AD[*:1])$min_AD" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
else
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.')" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
fi

tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.vcf.gz

echo $(gzip -d -c $out_dir/working_files/$prefix.subset.highquality.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Quality filtering complete'
