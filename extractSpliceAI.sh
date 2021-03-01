#!/bin/bash

while getopts "f:" opt; do
case $opt in
    f) input="$OPTARG";; # Input VCF file
esac
done

cp annotations/introme_annotate.vcf introme_annotate.spliceai.vcf

while read line; do
    if [[ ${line:0:1} != "#" ]]; then
        spliceai_scores=$(echo "$line" | grep -c "SpliceAI=")
        if (( "$spliceai_scores" > 0 )); then
            spliceai=$(echo "$line" | grep -o "SpliceAI=.*" | cut -f1)
            commas=$(echo "$spliceai" | grep -c ",")

            AG=()
            AL=()
            DG=()
            DL=()
            AG_dp=()
            AL_dp=()
            DG_dp=()
            DL_dp=()

            entries=($(echo "$spliceai" | tr ',' '\n'))
            for (( i=0; i<=${#commas[@]}; i++)); do
                AG+=($(echo ${entries[$i]} | cut -f3 -d '|'))
                AL+=($(echo ${entries[$i]} | cut -f4 -d '|'))
                DG+=($(echo ${entries[$i]} | cut -f5 -d '|'))
                DL+=($(echo ${entries[$i]} | cut -f6 -d '|'))
                AG_dp+=($(echo ${entries[$i]} | cut -f7 -d '|'))
                AL_dp+=($(echo ${entries[$i]} | cut -f8 -d '|'))
                DG_dp+=($(echo ${entries[$i]} | cut -f9 -d '|'))
                DL_dp+=($(echo ${entries[$i]} | cut -f10 -d '|'))
            done
            AG_max=$(echo ${AG[@]} | tr ' ' '\n' | sort -rn | head -1)
            AL_max=$(echo ${AL[@]} | tr ' ' '\n' | sort -rn | head -1)
            DG_max=$(echo ${DG[@]} | tr ' ' '\n' | sort -rn | head -1)
            DL_max=$(echo ${DL[@]} | tr ' ' '\n' | sort -rn | head -1)
            AG_dp_merge=$(echo ${AG_dp[@]} | tr ' ' ',')
            AL_dp_merge=$(echo ${AL_dp[@]} | tr ' ' ',')
            DG_dp_merge=$(echo ${DG_dp[@]} | tr ' ' ',')
            DL_dp_merge=$(echo ${DL_dp[@]} | tr ' ' ',')

            scores=$(echo "DS_AG=$AG_max;DS_AL=$AL_max;DS_DG=$DG_max;DS_DL=$DL_max;DP_AG=$AG_dp_merge;DP_AL=$AL_dp_merge;DP_DG=$DG_dp_merge;DP_DL=$DL_dp_merge")
            echo "$line" | awk -v s="$scores" '{$8=s; print }' OFS='\t' | cut -f 1-8 | grep -v "DS_AG=.;DS_AL=.;DS_DG=.;DS_DL=.;DP_AG=.;DP_AL=.;DP_DG=.;DP_DL=." >> introme_annotate.spliceai.vcf
        fi
    fi
done < $input

bcftools sort introme_annotate.spliceai.vcf -o introme_annotate.spliceai.vcf
uniq introme_annotate.spliceai.vcf | bgzip > introme_annotate.spliceai.vcf.gz
rm introme_annotate.spliceai.vcf
tabix -f introme_annotate.spliceai.vcf.gz
