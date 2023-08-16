#!/bin/bash

out_dir='./ese_test'
prefix='ese'
max_AF='0.01' # Minimum gnomAD allele frequency in the population with the highest allele frequency in gnomAD
reference_genome='AG_check/files/hg38_.fa'
# input="AG_check/files/pedcbioportal_long_44k.vcf"
input=$1

rm -rf $out_dir/working_files/
mkdir -p $out_dir/working_files/

##################################
# STEP 4: Hard filtering on the values of annotations added in the previous step

echo $(date +%x_%r) 'Beginning filtering by annotation values'

bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(gnomAD_PM_AF<=$max_AF || gnomAD_PM_AF='.')" $input | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz
tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz

variant_count=$(bcftools view -H $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | wc -l | tr -d ' ')
echo $(date +%x_%r) 'Frequency filtering complete -' $variant_count 'variants remaining'

if [[ $variant_count == 0 ]]; then
    exit 1
fi
 
# Remove annotations which interfere with score extraction
# bcftools view -h $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
grep -v "^#" $input | awk '{$8="."; print }' OFS='\t' >> $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
bgzip -f $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
gunzip -k $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf

echo $(date +%x_%r) 'Running ESE Scoring'
python3 ESE_ESS_scoring.py $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf $out_dir/working_files/$prefix.subset.highquality.ESE.vcf $reference_genome
cat annotations/introme_annotate.vcf $out_dir/working_files/$prefix.subset.highquality.ESE.vcf | bgzip > introme_annotate.ESE.vcf

