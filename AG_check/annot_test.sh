#!/bin/bash
# Tests the intermediary vcf file in the shellscript and python version for differences.


basedir=$(dirname "$0")

red=$'\e[1;31m'
grn=$'\e[1;32m'
end=$'\e[0m'

# for i in {1..2}; do
    for f in $basedir/files/pedcbioportal_strand_*; do
        echo -n "starting test on $f w/ python..."
        time python3 $basedir/AG_check.py $f $basedir/files/hg38_.fa;
        echo -n "starting test on $f w/ bash..."
        time bash $basedir/../AG_check.sh -i $f -r $basedir/files/hg38_.fa;

        bashresult="$basedir/../introme_annotate.functions.vcf.gz.tbi"
        pyresult="$basedir/../introme_annotate.functions2.vcf.gz.tbi"
        
        echo
        
        diff $bashresult $pyresult > /dev/null 2>&1
        
        error=$?
        if [ $error -eq 0 ]; then
            printf "%s\n" "${grn}Output vcfs are the same file${end}"
        elif [ $error -eq 1 ]; then
            printf "%s\n" "${red}Output vcfs differ${end}"
        else
            printf "%s\n" "${red}There was something wrong with the diff command${end}"
        fi

        echo "========================================="
    done;
# done;
