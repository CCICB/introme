#!/usr/bin/env bash

while getopts "r:v:" opt; do
case $opt in
    r) reference_genome="$OPTARG";; # output prefix
    v) input_vcf="$OPTARG";; # output prefix
esac
done

cp annotation_header.vcf introme_annotate.functions.vcf

while read line; do
    chr=$(echo "$line" | cut -f 1)
    pos=$(echo "$line" | cut -f 2)
    ref=$(echo "$line" | cut -f 4)
    alt=$(echo "$line" | cut -f 5)
    strand=$(echo "$line" | cut -f8 | grep -wo "Gene_Strand=.*" | sed 's/Gene_Strand=//g' | cut -f1 -d';')
    
    # Assign variant type for future annotations
    if [[ "${#ref}" == 1 ]] && [[ "${#alt}" == 1 ]]; then
        variant_type="SNV"
    elif [[ "${#ref}" > 1 ]] && [[ "${#alt}" == 1 ]]; then
        variant_type="INDEL"
    elif [[ "${#ref}" == 1 ]] && [[ "${#alt}" > 1 ]]; then
        variant_type="INDEL"
    elif [[ "${#ref}" > 1 ]] && [[ "${#alt}" > 1 ]]; then
        variant_type="INSDEL"
    fi

    # Define reference and alternative sequence
    ref_seq=$(samtools faidx $reference_genome "$chr:$(($pos-1))-$(($pos+${#ref}))" | grep -v "^>" | tr -d '\n')
    alt_seq=${ref_seq:0:1}$alt${ref_seq:$((${#ref}+1)):1}

    pos_strand=$(echo $strand | grep -c "+")
    neg_strand=$(echo $strand | grep -c "-")
    AG_found=0
    AG_lost=0
    GT_found=0
    GT_lost=0

    append=""

    if (( "$pos_strand" > 0 )); then
        if [[ "$alt_seq" == *"AG"* && "$ref_seq" != *"AG"* ]]; then
            append=$append"ag_created=+;"
        elif [[ "$alt_seq" != *"AG"* && "$ref_seq" == *"AG"* ]]; then
            append=$append"ag_lost=+;"
        elif [[ "$alt_seq" == *"GT"* && "$ref_seq" != *"GT"* ]]; then
            append=$append"gt_created=+;"
        elif [[ "$alt_seq" != *"GT"* && "$ref_seq" == *"GT"* ]]; then
            append=$append"gt_lost=+;"
        fi
    fi

    if (( "$neg_strand" > 0 )); then
        ref_seq_rev=$(echo "$ref_seq" | tr "[ACGT]" "[TGCA]" | rev)
        alt_seq_rev=$(echo "$alt_seq" | tr "[ACGT]" "[TGCA]" | rev)
        if [[ "$alt_seq_rev" == *"AG"* && "$ref_seq_rev" != *"AG"* ]]; then
            append=$append"ag_created=-;"
        elif [[ "$alt_seq_rev" != *"AG"* && "$ref_seq_rev" == *"AG"* ]]; then
            append=$append"ag_lost=-;"
        elif [[ "$alt_seq_rev" == *"GT"* && "$ref_seq_rev" != *"GT"* ]]; then
            append=$append"gt_created=-;"
        elif [[ "$alt_seq_rev" != *"GT"* && "$ref_seq_rev" == *"GT"* ]]; then
            append=$append"gt_lost=-;"
        fi
    fi

    echo "$chr"$'\t'"$pos"$'\t'"."$'\t'"$ref"$'\t'"$alt"$'\t'"1"$'\t'"PASS"$'\t'"variant_type=$variant_type;$append" >> introme_annotate.functions.vcf
done < <(bcftools view -H $input_vcf)

bcftools sort introme_annotate.functions.vcf -o introme_annotate.functions.vcf
uniq introme_annotate.functions.vcf | bgzip > introme_annotate.functions.vcf.gz
rm introme_annotate.functions.vcf
tabix -f introme_annotate.functions.vcf.gz
