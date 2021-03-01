#!/usr/bin/env bash

while getopts "p:r:" opt; do
case $opt in
    p) prefix="$OPTARG";; # output prefix
    r) reference_genome="$OPTARG";; # output prefix
esac
done

cp annotations/introme_annotate.vcf introme_annotate.functions.vcf

bcftools view -H output/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | \
while read line; do
    chr=$(echo "$line" | cut -f 1)
    pos=$(echo "$line" | cut -f 2)
    id=$(echo "$line" | cut -f 3)
    ref=$(echo "$line" | cut -f 4)
    alt=$(echo "$line" | cut -f 5)
    strand=$(echo "$line" | cut -f8 | grep -wo "Gene_Strand=.*" | sed 's/Gene_Strand=//g' | cut -f1 -d';')

    if [[ "${#ref}" == 1 ]] && [[ "${#alt}" == 1 ]]; then
        variant_type="SNV"
    elif [[ "${#ref}" > 1 ]] && [[ "${#alt}" == 1 ]]; then
        variant_type="INDEL"
    elif [[ "${#ref}" == 1 ]] && [[ "${#alt}" > 1 ]]; then
        variant_type="INDEL"
    elif [[ "${#ref}" > 1 ]] && [[ "${#alt}" > 1 ]]; then
        variant_type="INSDEL"
    fi

    in_AG_region=$(echo "$line" | cut -f8 | grep -c -e "branchpoint" -e "acceptor_region" -e "acceptor_canonical")
    in_GT_region=$(echo "$line" | cut -f8 | grep -v -e "exonic" -e "canonical" -e "branchpoint" -e "acceptor" | grep -c -e "5_intronic" -e "donor_region")
    if [[ "$in_AG_region" == 1 || "$in_GT_region" == 1 ]]; then
        ref_seq=$(samtools faidx $reference_genome "$chr:$(($pos-1))-$(($pos+${#ref}))" | grep -v "^>")
        alt_seq=${ref_seq:0:1}$alt${ref_seq:$((${#ref}+1)):1}

        pos_strand=$(echo $strand | grep -c "+")
        neg_strand=$(echo $strand | grep -c "-")
        AG_found=0
        GT_found=0

        if (( "$pos_strand" > 0 )); then
            if [[ "$in_AG_region" == 1 && "$alt_seq" == *"AG"* ]]; then
                AG_found=1
            elif [[ "$in_GT_region" == 1 && "$alt_seq" == *"GT"* ]]; then
                GT_found=1
            fi
        fi

        if (( "$neg_strand" > 0 )); then
            alt_seq_rev=$(echo "$alt_seq" | tr "[ACGT]" "[TGCA]" | rev)
            if [[ "$in_AG_region" == 1 && "$alt_seq_rev" == *"AG"* ]]; then
                AG_found=1
            elif [[ "$in_GT_region" == 1 && "$alt_seq_rev" == *"GT"* ]]; then
                GT_found=1
            fi
        fi

        if [[ "$AG_found" == 1 ]]; then
            echo "$chr"$'\t'"$pos"$'\t'"$id"$'\t'"$ref"$'\t'"$alt"$'\t'"1"$'\t'"PASS"$'\t'"variant_type=$variant_type;ag_created=1" >> introme_annotate.functions.vcf
        elif [[ "$GT_found" == 1 ]]; then
            echo "$chr"$'\t'"$pos"$'\t'"$id"$'\t'"$ref"$'\t'"$alt"$'\t'"1"$'\t'"PASS"$'\t'"variant_type=$variant_type;gt_created=1" >> introme_annotate.functions.vcf
        else
            echo "$chr"$'\t'"$pos"$'\t'"$id"$'\t'"$ref"$'\t'"$alt"$'\t'"1"$'\t'"PASS"$'\t'"variant_type=$variant_type" >> introme_annotate.functions.vcf
        fi
    else
        echo "$chr"$'\t'"$pos"$'\t'"$id"$'\t'"$ref"$'\t'"$alt"$'\t'"1"$'\t'"PASS"$'\t'"variant_type=$variant_type" >> introme_annotate.functions.vcf
    fi
done

bcftools sort introme_annotate.functions.vcf -o introme_annotate.functions.vcf
uniq introme_annotate.functions.vcf | bgzip > introme_annotate.functions.vcf.gz
rm introme_annotate.functions.vcf
tabix -f introme_annotate.functions.vcf.gz
