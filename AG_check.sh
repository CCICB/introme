#!/usr/bin/env bash

while getopts "i:r:" opt; do
case $opt in
    i) input_file="$OPTARG";; # Input file
    r) reference_genome="$OPTARG";; # Reference Genome
esac
done

docker run -v $(pwd):/app ag_check \
    python3 AG_check/AG_check.py \
    $input_file \
    $reference_genome

bcftools sort introme_annotate.functions.vcf | uniq | bgzip > introme_annotate.functions.vcf.gz
rm introme_annotate.functions.vcf
tabix -f introme_annotate.functions.vcf.gz
