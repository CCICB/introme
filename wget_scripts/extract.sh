#!/bin/bash

# Define inputs for score extraction portion of introme
while getopts "a:b:c:p:" opt; do
    case $opt in
        a) spliceai_vcf="$OPTARG";; # SpliceAI output VCF
        b) spliceogen_txt="$OPTARG";; # Spliceogen output txt
        c) mmsplice_vcf="$OPTARG";; # MMSplice output vcf
        p) prefix="$OPTARG";; # Output file prefix
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list


echo $(date +%x_%r) 'Beginning SpliceAI score extraction'

cp annotation_header.vcf introme_annotate.spliceai.vcf
bcftools view -h $spliceai_vcf > $prefix.spliceai.rmcommas.vcf
bcftools view -H $spliceai_vcf | awk '$8 !~ /,/' >> $prefix.spliceai.rmcommas.vcf

# If there are multiple scores per line, keep the highest values and format correctly
while read line; do
    start=$(echo "$line" | cut -f1-7)
    end=$(echo "$line" | cut -f9-11)

    scores_gene1=$(echo "$line" | cut -f8 | cut -f3-6 -d'|')
    scores_gene2=$(echo "$line" | cut -f8 | cut -f12-15 -d'|')

    if [[ $scores_gene1 == $scores_gene2 ]]; then
        gene_concat=$(echo "$line" | awk -F "|" '{ print $2 "," $11 }')
        spliceai_scores=$(echo "$line" | cut -f8 | awk -F'|' -v t="$gene_concat" '{$2=t; print }' OFS='|' | cut -f1,2 -d',')
        echo "$start"$'\t'"$spliceai_scores"$'\t'"$end" >> $prefix.spliceai.rmcommas.vcf
    elif [[ $(echo "$line" | grep "0.00|0.00|0.00|0.00") ]]; then
        echo "$(echo "$start" | tr '\n' '\t' && echo "$line" | cut -f8 | tr ',' '\n' | grep -v "0.00|0.00|0.00|0.00" | sed 's/|/}/' | sed 's/.*}/SpliceAI=.|/' | tr '\n' '\t' && echo "$end")" >> $prefix.spliceai.rmcommas.vcf
    else
        max_gene1=$(echo "$scores_gene1" | tr '|' '\n' | sort -rn | head -1)
        max_gene2=$(echo "$scores_gene2" | tr '|' '\n' | sort -rn | head -1)
        if (( $(echo "$max_gene1 >= $max_gene2" | bc -l) )); then
            spliceai_scores=$(echo "$line" | cut -f8 | cut -f1 -d',')
        elif (( $(echo "$max_gene1 < $max_gene2" | bc -l) )); then
            spliceai_scores=$(echo "$line" | cut -f8 | cut -f2 -d',' | sed 's/|/}/' | sed 's/.*}/SpliceAI=.|/')
        fi
        echo "$start"$'\t'"$spliceai_scores"$'\t'"$end" >> $prefix.spliceai.rmcommas.vcf
    fi
done < <(bcftools view -H $spliceai_vcf | awk '$8 ~ /,/')

# Format SpliceAI output for vcfanno
bcftools view -H $prefix.spliceai.rmcommas.vcf | cut -f1-8 | grep "SpliceAI" | sed 's/|/gene=/' | sed 's/SpliceAI=*.*gene/gene/' | sed 's/|/;DS_AG=/' | sed 's/|/;DS_AL=/' | sed 's/|/;DS_DG=/' | sed 's/|/;DS_DL=/' | sed 's/|/;DP_AG=/' | sed 's/|/;DP_AL=/' | sed 's/|/;DP_DG=/' | sed 's/|/;DP_DL=/' >> introme_annotate.spliceai.vcf
bcftools sort introme_annotate.spliceai.vcf -o introme_annotate.spliceai.sort.vcf
uniq introme_annotate.spliceai.sort.vcf | bgzip > introme_annotate.spliceai.vcf.gz
tabix -f introme_annotate.spliceai.vcf.gz



echo $(date +%x_%r) 'Beginning Spliceogen score extraction'

# Format Spliceogen txt file to a vcf file and relabel for annotation uses
cp annotation_header.vcf introme_annotate.spliceogen.vcf
grep -v "^#" $spliceogen_txt | \
  cut -f1,2,4,5,8-11,16-23 | sed -e $'s/\t/\t.\t/2' -e $'s/\t/\t.\tPASS\t/5' \
  -e $'s/\t/\tmesDonRef=/7' -e $'s/\t/;mesDonAlt=/8' -e $'s/\t/;mesAccRef=/8' \
  -e $'s/\t/;mesAccAlt=/8' -e $'s/\t/;ESEmaxRef=/8' -e $'s/\t/;ESEmaxAlt=/8' \
  -e $'s/\t/;ESEminRef=/8' -e $'s/\t/;ESEminAlt=/8' -e $'s/\t/;DonGainP=/8' \
  -e $'s/\t/;AccGainP=/8' -e $'s/\t/;DonLossP=/8' -e $'s/\t/;AccLossP=/8' \
    >> introme_annotate.spliceogen.vcf
bcftools sort introme_annotate.spliceogen.vcf -o introme_annotate.spliceogen.sort.vcf
uniq introme_annotate.spliceogen.sort.vcf | bgzip > introme_annotate.spliceogen.vcf.gz
tabix introme_annotate.spliceogen.vcf.gz



echo $(date +%x_%r) 'Beginning MMSplice score extraction'

# Extract MMSplice Scores
cp annotation_header.vcf introme_annotate.mmsplice.vcf
bcftools view -H $mmsplice_vcf | cut -f1-8 | grep "CSQ" | \
  sed -e 's/CSQ=/alt_acceptor=/' -e 's/|/;alt_acceptor_intron=/' -e 's/|/;alt_donor=/' \
  -e 's/|/;alt_donor_intron=/' -e 's/|/;alt_exon=/' -e 's/|/;delta_logit_PSI=/' \
  -e 's/|/;pathogenicity=/' -e 's/|/;ref_acceptor=/' -e 's/|/;ref_acceptor_intron=/' \
  -e 's/|/;ref_donor=/' -e 's/|/;ref_donor_intron=/' -e 's/|/;ref_exon=/' \
        >> introme_annotate.mmsplice.vcf
bcftools sort introme_annotate.mmsplice.vcf -o introme_annotate.mmsplice.sort.vcf
uniq introme_annotate.mmsplice.sort.vcf | bgzip > introme_annotate.mmsplice.vcf.gz
tabix introme_annotate.mmsplice.vcf.gz


echo $(date +%x_%r) 'Finished score extraction'
