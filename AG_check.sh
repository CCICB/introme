#!/usr/bin/env bash

while getopts "i:r:" opt; do
case $opt in
    i) input_file="$OPTARG";; # Input file
    r) reference_genome="$OPTARG";; # Reference Genome
esac
done

cp annotations/introme_annotate.vcf introme_annotate.functions.vcf

while read line; do
    chr=$(echo "$line" | cut -f 1)
    pos=$(echo "$line" | cut -f 2)
    ref=$(echo "$line" | cut -f 4)
    alt=$(echo "$line" | cut -f 5)
    strand=$(echo "$line" | cut -f8 | grep -wo "strand=.*" | sed 's/strand=//g' | cut -f1 -d';')

    if [[ "${#ref}" == 1 ]] && [[ "${#alt}" == 1 ]]; then
        variant_type="SNV"
    elif [[ "${#ref}" > 1 ]] && [[ "${#alt}" == 1 ]]; then
        variant_type="INDEL"
    elif [[ "${#ref}" == 1 ]] && [[ "${#alt}" > 1 ]]; then
        variant_type="INDEL"
    elif [[ "${#ref}" > 1 ]] && [[ "${#alt}" > 1 ]]; then
        variant_type="INSDEL"
    fi

    # 1-based
    ref_seq=$(samtools faidx $reference_genome "$chr:$(($pos-1))-$(($pos+${#ref}))" | grep -v "^>" | tr -d '\n')
    alt_seq=${ref_seq:0:1}$alt${ref_seq:$((${#ref}+1)):1}

    # echo "$pos, $ref_seq, $alt_seq, $ref, $alt"

    pos_strand=$(echo $strand | grep -c "+")
    neg_strand=$(echo $strand | grep -c "-")

    append=""

    if (( "$pos_strand" > 0 )); then
        if [[ "$alt_seq" == *"AG"* && "$ref_seq" != *"AG"* ]]; then
            append+="ag_created=+;"
        fi
        if [[ "$alt_seq" != *"AG"* && "$ref_seq" == *"AG"* ]]; then
            append+="ag_lost=+;"
        fi
        if [[ "$alt_seq" == *"GT"* && "$ref_seq" != *"GT"* ]]; then
            append+="gt_created=+;"
        fi
        if [[ "$alt_seq" != *"GT"* && "$ref_seq" == *"GT"* ]]; then
            append+="gt_lost=+;"
        fi
    fi

    if (( "$neg_strand" > 0 )); then
        ref_seq_rev=$(echo "$ref_seq" | tr "[ACGT]" "[TGCA]" | rev)
        alt_seq_rev=$(echo "$alt_seq" | tr "[ACGT]" "[TGCA]" | rev)
        if [[ "$alt_seq_rev" == *"AG"* && "$ref_seq_rev" != *"AG"* ]]; then
            append+="ag_created=-;"
        fi
        if [[ "$alt_seq_rev" != *"AG"* && "$ref_seq_rev" == *"AG"* ]]; then
            append+="ag_lost=-;"
        fi
        if [[ "$alt_seq_rev" == *"GT"* && "$ref_seq_rev" != *"GT"* ]]; then
            append+="gt_created=-;"
        fi
        if [[ "$alt_seq_rev" != *"GT"* && "$ref_seq_rev" == *"GT"* ]]; then
            append+="gt_lost=-;"
        fi
    fi

    echo "$chr"$'\t'"$pos"$'\t'"."$'\t'"$ref"$'\t'"$alt"$'\t'"."$'\t'"PASS"$'\t'"variant_type=$variant_type;$append" >> introme_annotate.functions.vcf
done < <(bcftools view -H $input_file)

# bcftools sort introme_annotate.functions.vcf | uniq | bgzip > introme_annotate.functions.vcf.gz
# rm introme_annotate.functions.vcf
# tabix -f introme_annotate.functions.vcf.gz
