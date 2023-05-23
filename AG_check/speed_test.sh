#!/bin/bash
basedir=$(dirname "$0")

for f in $basedir/files/pedcbioportal_long_*; do
    for i in {1..2}; do
        echo "==========================="
        echo -n "starting test on $f w/ python..."
        time python3 $basedir/AG_check.py $f $basedir/files/hg38_.fa 1>/dev/null;

        echo
    done;
done;
