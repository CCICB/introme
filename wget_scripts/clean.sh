#!/bin/bash

# Define inputs for cleaning portion of introme
while getopts "p:v:" opt; do
    case $opt in
        p) prefix="$OPTARG";;
        v) annotated_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

while read line; do
    if [[ ${line:0:1} == "#" ]]; then
        echo "$line" > $prefix.introme.tsv
    else
        exonic_count=$(echo "$line" | cut -f 10 | grep -c "exonic")
        intronic_count=$(echo "$line" | cut -f 10 | grep -c "intronic")
        if (( $exonic_count >= 1 )) && (( $intronic_count == 0 )); then
            gene_location="exonic"
        elif (( $intronic_count >= 1 )) && (( $exonic_count == 0 )); then
            gene_location="intronic"
        elif (( $intronic_count >= 1 )) && (( $exonic_count >= 1)); then
            gene_location="exonic,intronic"
        else
            gene_location="."
        fi

        variant_type=$(echo "$line" | cut -f 7 | cut -f1 -d',' )

        regions=$(echo "$line" | cut -f 12)
        if [[ $regions == "." || $regions != *","* ]]; then
            region_replace=$regions
        elif [[ $regions == *"donor_canonical"* ]]; then
            region_replace="donor_canonical"
        elif [[ $regions == *"acceptor_canonical"* ]]; then
            region_replace="acceptor_canonical"
        elif [[ $regions == *"donor_region"* ]]; then
            region_replace="donor_region"
        elif [[ $regions == *"acceptor_region"* ]]; then
            region_replace="acceptor_region"
        elif [[ $regions == *"donor_exonic"* ]]; then
            region_replace="donor_exonic"
        elif [[ $regions == *"acceptor_exonic"* ]]; then
            region_replace="acceptor_exonic"
        elif [[ $regions == *"branchpoint"* ]]; then
            region_replace="branchpoint"
        elif [[ $regions == *"5_intronic"* ]]; then
            region_replace="5_intronic"
        fi

        echo "$line" | awk -v v="$variant_type" -v g="$gene_location" -v r="$region_replace" '{$7=v; $10=g; $12=r; print }' OFS='\t' >> $prefix.introme.tsv
    fi
done < $annotated_vcf
