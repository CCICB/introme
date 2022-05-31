#!/bin/bash

# Define inputs for cleaning portion of introme
while getopts "p:v:" opt; do
    case $opt in
        p) prefix="$OPTARG";;
        v) annotated_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

if [[ ! -z $(bcftools view $annotated_VCF | grep "CSQ") ]]; then
    bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/Variant_Type\t%INFO/Gene_Symbol\t%INFO/Gene_Strand\t%INFO/Gene_Location\t%INFO/Intron_Type\t%INFO/Gene_Regions\t%INFO/CADD_Phred\t%INFO/MGRB_AF\t%INFO/gnomAD_PM_AF\t%INFO/CSQ\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/mesAccRef\t%INFO/mesAccAlt\t%INFO/mesDonRef\t%INFO/mesDonAlt\t%INFO/AG_Created\t%INFO/AG_Lost\t%INFO/GT_Created\t%INFO/GT_Lost\t%INFO/Branchpointer_Prob\t%INFO/Branchpointer_U2_Binding_Energy\t%INFO/Branchpointer_options\t%INFO/Branchpointer_max_U2_Binding_Energy\t%INFO/Branchpointer_max_Prob\t%INFO/SpliceAI_Gene\t%INFO/SpliceAI_Acceptor_Gain\t%INFO/SpliceAI_Acceptor_Loss\t%INFO/SpliceAI_Donor_Gain\t%INFO/SpliceAI_Donor_Loss\t%INFO/SpliceAI_DP_Acceptor_Gain\t%INFO/SpliceAI_DP_Acceptor_Loss\t%INFO/SpliceAI_DP_Donor_Gain\t%INFO/SpliceAI_DP_Donor_Loss\t%INFO/AccGainP\t%INFO/AccLossP\t%INFO/DonGainP\t%INFO/DonLossP\t%INFO/MMSplice_alt_acceptor\t%INFO/MMSplice_alt_acceptor_intron\t%INFO/MMSplice_alt_donor\t%INFO/MMSplice_alt_donor_intron\t%INFO/MMSplice_alt_exon\t%INFO/MMSplice_delta_logit_PSI\t%INFO/MMSplice_pathogenicity\t%INFO/MMSplice_ref_acceptor\t%INFO/MMSplice_ref_acceptor_intron\t%INFO/MMSplice_ref_donor\t%INFO/MMSplice_ref_donor_intron\t%INFO/MMSplice_ref_exon\t%INFO/ESEmaxRef\t%INFO/ESEmaxAlt\t%INFO/ESEminRef\t%INFO/ESEminAlt[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $annotated_VCF > $prefix.annotated.tsv
else
    bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/Variant_Type\t%INFO/Gene_Symbol\t%INFO/Gene_Strand\t%INFO/Gene_Location\t%INFO/Intron_Type\t%INFO/Gene_Regions\t%INFO/CADD_Phred\t%INFO/MGRB_AF\t%INFO/gnomAD_PM_AF\t.\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/mesAccRef\t%INFO/mesAccAlt\t%INFO/mesDonRef\t%INFO/mesDonAlt\t%INFO/AG_Created\t%INFO/AG_Lost\t%INFO/GT_Created\t%INFO/GT_Lost\t%INFO/Branchpointer_Prob\t%INFO/Branchpointer_U2_Binding_Energy\t%INFO/Branchpointer_options\t%INFO/Branchpointer_max_U2_Binding_Energy\t%INFO/Branchpointer_max_Prob\t%INFO/SpliceAI_Gene\t%INFO/SpliceAI_Acceptor_Gain\t%INFO/SpliceAI_Acceptor_Loss\t%INFO/SpliceAI_Donor_Gain\t%INFO/SpliceAI_Donor_Loss\t%INFO/SpliceAI_DP_Acceptor_Gain\t%INFO/SpliceAI_DP_Acceptor_Loss\t%INFO/SpliceAI_DP_Donor_Gain\t%INFO/SpliceAI_DP_Donor_Loss\t%INFO/AccGainP\t%INFO/AccLossP\t%INFO/DonGainP\t%INFO/DonLossP\t%INFO/MMSplice_alt_acceptor\t%INFO/MMSplice_alt_acceptor_intron\t%INFO/MMSplice_alt_donor\t%INFO/MMSplice_alt_donor_intron\t%INFO/MMSplice_alt_exon\t%INFO/MMSplice_delta_logit_PSI\t%INFO/MMSplice_pathogenicity\t%INFO/MMSplice_ref_acceptor\t%INFO/MMSplice_ref_acceptor_intron\t%INFO/MMSplice_ref_donor\t%INFO/MMSplice_ref_donor_intron\t%INFO/MMSplice_ref_exon\t%INFO/ESEmaxRef\t%INFO/ESEmaxAlt\t%INFO/ESEminRef\t%INFO/ESEminAlt[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $annotated_VCF > $prefix.annotated.tsv
fi

# Remove the '[<number]' column numbers added by bcftools query to column names
grep '^#' $prefix.annotated.tsv | sed "s/\[[0-9]*\]//g" > $prefix.annotated.cleaned.tsv
grep -v '^#' $prefix.annotated.tsv >> $prefix.annotated.cleaned.tsv

# Sort by chromosome and coordinate
sort -k1,1n -k2,2n $prefix.annotated.cleaned.tsv > $prefix.annotated.tsv

while read line; do
    if [[ ${line:0:1} == "#" ]]; then
        echo "$line" > $prefix.annotated.cleaned.tsv
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

        echo "$line" | awk -v v="$variant_type" -v g="$gene_location" -v r="$region_replace" '{$7=v; $10=g; $12=r; print }' OFS='\t' >> $prefix.annotated.cleaned.tsv
    fi
done < $prefix.annotated.tsv
