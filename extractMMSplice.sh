#!/bin/bash

while getopts "f:p:" opt; do
case $opt in
    f) file="$OPTARG";; # Input file of scores
    p) prefix="$OPTARG";; # output prefix
esac
done

sed -n '/##fileformat=VCF/,$p' $file > output/working_files/$prefix.mmsplice.vcf
cp annotations/introme_annotate.vcf introme_annotate.mmsplice.vcf

# bcftools view -h output/working_files/$prefix.mmsplice.vcf > output/working_files/$prefix.mmsplice.fix.vcf
# bcftools view -H output/working_files/$prefix.mmsplice.vcf | \
# while read line;
# do
#     mmsplice=$(echo "$line" | cut -f8 | grep -wo "CSQ=.*")
#     echo "$line" | awk -v t="$mmsplice" '{$8=t; print }' OFS='\t' >> output/working_files/$prefix.mmsplice.fix.vcf
# done
#
# mv output/working_files/$prefix.mmsplice.fix.vcf output/working_files/$prefix.mmsplice.vcf


bcftools view -H output/working_files/$prefix.mmsplice.vcf | \
while read line;
do
    mmsplice=$(echo "$line" | cut -f8 | grep -o "CSQ=.*")
    commas=$(echo $mmsplice | grep "MM" | tr -dc "," | wc -c )
    #commas=$(echo $line | cut -f8 | grep "MM" | tr -dc "," | wc -c )
    entries=$(($commas+2))
    logit=()
    score_comparison=""
    highest_score="."

    for (( i = 1; i < "$entries"; i++ ));
    do
        logit+=($(echo $mmsplice | cut -d'"' -f2 | cut -d',' -f"$i" | cut -d'|' -f29))
        logit+=($(echo $line | cut -d'"' -f2 | cut -d',' -f"$i" | cut -d'|' -f29))

        current_score=$(echo $line | cut -d'"' -f2 | cut -d',' -f"$i" | cut -d'|' -f29)
        abs_current_score=$(echo "$current_score" | tr -d "-")

        if [[ "$abs_current_score">"$score_comparison" ]]; then
            score_comparison=$abs_current_score
            highest_score=$current_score
        fi
    done

    if [[ $highest_score != "." ]]; then
        rounded=$(printf '%.*f\n' 2 $highest_score)
        mmsplice=$(echo "mmsplice_dl="$rounded)
        echo "$line" | awk -v t="$mmsplice" '{$8=t; print }' OFS='\t' | cut -f 1-8 >> introme_annotate.mmsplice.vcf
    fi
done

bcftools sort introme_annotate.mmsplice.vcf -o introme_annotate.mmsplice.vcf
uniq introme_annotate.mmsplice.vcf | bgzip > introme_annotate.mmsplice.vcf.gz
rm introme_annotate.mmsplice.vcf
tabix -f introme_annotate.mmsplice.vcf.gz
