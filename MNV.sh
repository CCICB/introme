#!/usr/bin/env bash


while getopts "a:f:p:r:" opt; do
    case $opt in
        a) genome="$OPTARG";; # Genome build
        f) file="$OPTARG";; # Original VCF
        p) prefix="$OPTARG";;
        r) reference_genome="$OPTARG";;
    esac
done

bcftools filter -i"TYPE!='snp' && TYPE!='indel'" $file | grep -v "^#" > output/working_files/$prefix.mnv.vars.vcf

if [[ $genome == "hg19" ]]; then
    assembly="grch37"
elif [[ $genome == "hg38" ]]; then
    assembly="grch38"
fi

while read line; do
    bcftools view -h $file > output/working_files/$prefix.mnv.vcf
    chr=$(echo "$line" | cut -f1)
    pos=$(echo "$line" | cut -f2)
    ref=$(echo "$line" | cut -f4)
    alt=$(echo "$line" | cut -f5)
    end=$(echo "$line" | cut -f6-11)

    echo "$chr"$'\t'"$pos"$'\t'"del"$'\t'"$ref"$'\t'"${ref:0:1}"$'\t'"$end" >> output/working_files/$prefix.mnv.vcf
    echo "$chr"$'\t'"$pos"$'\t'"ins1"$'\t'"${ref:0:1}"$'\t'"$alt"$'\t'"$end" >> output/working_files/$prefix.mnv.vcf
    echo "$chr"$'\t'"$(($pos+${#ref}-1))"$'\t'"ins2"$'\t'"${ref:$((${#ref}-1)):${#ref}}"$'\t'"${ref:$((${#ref}-1)):${#ref}}""${alt:1:${#alt}}"$'\t'"$end" >> output/working_files/$prefix.mnv.vcf

    spliceAI=$(spliceai -I output/working_files/$prefix.mnv.vcf -A $assembly -R $reference_genome -D 1000 | grep -v "^#" | cut -f 8)

    AG=$(echo "$spliceAI" | cut -f3 -d'|' | sort -nr | head -1)
    AL=$(echo "$spliceAI" | cut -f4 -d'|' | sort -nr | head -1)
    DG=$(echo "$spliceAI" | cut -f5 -d'|' | sort -nr | head -1)
    DL=$(echo "$spliceAI" | cut -f6 -d'|' | sort -nr | head -1)

    start=$(echo "$line" | cut -f1-7)
    echo "$start"$'\t'"DS_AG=$AG;DS_AL=$AL;DS_DG=$DG;DS_DL=$DL" >> introme_annotate.spliceai.vcf

done < output/working_files/$prefix.mnv.vars.vcf

rm output/working_files/$prefix.mnv.vcf

bcftools sort introme_annotate.spliceai.vcf -Oz -o introme_annotate.spliceai.vcf.gz
rm introme_annotate.spliceai.vcf
tabix -f introme_annotate.spliceai.vcf.gz
