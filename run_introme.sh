#!/bin/bash

# Introme Bash Script
#
# Developers:
# Patricia Sullivan (psullivan@ccia.org.au)
# Velimir Gayevskiy (vel@vel.nz)
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)


# Version=0.4.1
##################################
# DEFINE VARIABLES

# Directories
out_dir='./output'

# Trigger values
help=0
no_csq=0
no_pop_filter=0
no_qual_filter=0
mode=full

# Hard filter cutoffs for each annotation, filtering thresholds will depend on your application, we urge you to carefully consider these
MGRB_AF='<=0.01' # Minimum MGRB allele frequency (healthy Australian population)
gnomad_popmax_AF='<=0.01' # Minimum gnomAD allele frequency in the population with the highest allele frequency in gnomAD
CADD_Phred='>=0' # Minimum CADD score (phred-scaled)
min_QUAL='>=200' # The QUAL VCF field
min_DP='>=20' # The sample with the highest depth must have a equal or larger depth than this value
min_AD='>=5' # The sample with the highest number of alternate reads must have a equal or larger number than this value

##################################
# INTROME ARGUMENTS

while getopts "b:cfg:hm:p:qr:v::" opt; do
    case $opt in
        b) input_BED="$OPTARG";; # Input BED file (i.e. regions of interest)
        c) no_csq=1;; # trigger no csq command
        f) no_pop_filter=1;; # trigger no population filter command
        g) gene_list+="$OPTARG";; # Input list of genes to score
        h) help=1;; # trigger help command
        m) mode="$OPTARG";; # Mode in which introme is running - add "fast" for precomputed scores only
        p) prefix="$OPTARG";; # Output file prefix
        q) no_qual_filter=1;; # trigger no quality filter command
        r) reference_genome="$OPTARG";; # Path to the reference genome used for mapping
        v) input_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

# Run help command
if [[ $help == 1 ]]; then
    echo "Program: Introme (A program for variant calling and manipulating VCFs and BCFs)"
    echo "Version: 0.4.1"$'\n'

    echo "Usage: ./run_introme.sh [options] -b <subset.bed.gz> -p <prefix> -r <reference_genome.fa> -v <variants.vcf.gz>"$'\n'

    echo "Required Commands:"
    echo "-b <subset.bed.gz>        Input BED file of the regions of interest"
    echo "-p <prefix>               Output file prefix"
    echo "-r <reference_genome.fa>  The reference genome used for mapping"
    echo "-v <variants.vcf.gz>      Input variants in VCF format (must be gzipped)"$'\n'

    echo "Options:"
    echo "-c                        Use if the input VCF file has not been run through VEP"
    echo "-f                        Use to remove allele frequency filtering (keep common variants)"
    echo "-g                        Genes of interest"
    echo "-h                        Print Introme usage instructions"
    echo "-m fast                   Run Introme fast mode (will only compute SNVs)"
    echo "-q                        Use to remove variant quality filtering"$'\n'
    exit 1
fi

# Make sure each required argument has a non-empty value
if [ -z $input_BED ]; then
	echo "No restriction BED file supplied."
	exit 1
elif [ -z $input_VCF ]; then
	echo "No path to an input VCF file has been supplied."
	exit 1
elif [ -z $prefix ]; then
	echo "No output prefix has been supplied."
	exit 1
elif [ -z $reference_genome ]; then
	echo "No reference genome file has been supplied."
	exit 1
fi


##################################
# STEP 1: subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)

echo $(date +%x_%r) 'Beginning subsetting to genomic regions of interest'
echo $(gzip -d -c $input_VCF | grep -v '^#' | wc -l) 'variants prior to subsetting'

bcftools sort $input_VCF -Oz -o $out_dir/working_files/$prefix.sorted.vcf.gz &>/dev/null # Ensures the file is sorted correctly prior to subsetting

if [ ! -z $gene_list ]; then
    for (( gene_entry = 0; gene_entry < "${#gene_list[@]}"; gene_entry++ ));
    do
        echo $(date +%x_%r) 'Gene list provided - Beginning subsetting to genes of interest'
        if (( $(echo "${gene_list[$gene_entry]}" | grep -c "\.txt") )); then
            grep -wf "${gene_list[$gene_entry]}" $input_BED >> subsetting/$prefix.genelist.bed
        else
            grep -w "${gene_list[$gene_entry]}" $input_BED >> subsetting/$prefix.genelist.bed
        fi
        bedtools intersect -header -u -a $out_dir/working_files/$prefix.sorted.vcf.gz -b subsetting/$prefix.genelist.bed | bgzip > $out_dir/working_files/$prefix.subset.vcf.gz # -u for unique record in VCF, otherwise multiple variants are output for overlapping introns
    done
elif [ ! -z $zero_categories ]; then
    echo $(date +%x_%r) 'ZERO category provided - Beginning subsetting to '$zero_categories' genes'
    bedtools intersect -header -u -a $out_dir/working_files/$prefix.sorted.vcf.gz -b ~/Projects/Landscape/Reported_BEDs/$zero_categories.bed | bgzip > $out_dir/working_files/$prefix.subset.vcf.gz
else
    echo $(date +%x_%r) 'No gene list provided - Beginning subsetting to genomic regions of interest'
    bedtools intersect -header -u -a $out_dir/working_files/$prefix.sorted.vcf.gz -b $input_BED | bgzip > $out_dir/working_files/$prefix.subset.vcf.gz # -u for unique record in VCF
fi

tabix -p vcf $out_dir/working_files/$prefix.subset.vcf.gz

echo $(gzip -d -c $out_dir/working_files/$prefix.subset.vcf.gz | grep -v '^#' | wc -l) 'variants after subsetting'
echo $(date +%x_%r) 'Subsetting complete'


# STEP 3: Hard filtering on variant quality (this is here to reduce the number of variants going into the CPU-costly annotation step below)

echo $(date +%x_%r) 'Beginning quality filtering'

if [ $no_qual_filter == 0 ]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.') && MAX(FORMAT/DP[*])$min_DP && MAX(FORMAT/AD[*:1])$min_AD" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
else
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.')" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
fi

tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.vcf.gz

echo $(gzip -d -c $out_dir/working_files/$prefix.subset.highquality.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Quality filtering complete'

##################################
# STEP 4: annotate the subsetted VCF with useful information, to be used for filtering downstream


export IRELATE_MAX_GAP=1000 # Note: this is set to speed up annotation when .csi (as opposed to .tbi) files are present via https://github.com/brentp/vcfanno/releases/

# Run the appropriate annotation command based off the mode supplied
if [ $mode == "full" ]; then
    echo $(date +%x_%r) 'Beginning annotation'
    vcfanno -p $(getconf _NPROCESSORS_ONLN) conf.toml $out_dir/working_files/$prefix.subset.highquality.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
elif [ $mode == "fast" ]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"TYPE='snp'" $out_dir/working_files/$prefix.subset.highquality.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.SNPs.vcf.gz # Filter file for SNPs
    tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.SNPs.vcf.gz
    echo $(gzip -d -c $out_dir/working_files/$prefix.subset.highquality.SNPs.vcf.gz | grep -v '^#' | wc -l) 'variants after SNP filtering for fast mode'

    echo $(date +%x_%r) 'Beginning annotation'
    vcfanno -p $(getconf _NPROCESSORS_ONLN) conf_fast.toml $out_dir/working_files/$prefix.subset.highquality.SNPs.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
fi

tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz

echo $(date +%x_%r) 'Annotation complete'

##################################
# STEP 5: Hard filtering on the values of annotations added in the previous step

echo $(date +%x_%r) 'Beginning filtering by annotation values'

if [ $no_pop_filter == 0 ]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(MGRB_AF$MGRB_AF || MGRB_AF='.') && (gnomAD_PM_AF$gnomad_popmax_AF || gnomAD_PM_AF='.')" $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz
else
    cp $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz
fi
tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz


echo $(gzip -d -c $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | grep -v '^#' | wc -l) 'variants after filtering'
echo $(date +%x_%r) 'Frequency filtering complete'

##################################
# STEP 6: Run MMSplice and Splice AI - Skipped if fast mode

if [ $mode == "full" ]; then
    gunzip -k $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz # MMSplice and SpliceAI needs unzipped input files

    # Remove annotations which interfere with score extraction
    bcftools view -h $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
    grep -v "^#" $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf | awk '{$8="."; print }' OFS='\t' >> $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
    mv $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf

    # Run MMSplice using docker
    echo $(date +%x_%r) 'Starting to run MMSplice'
    cat $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf | sudo docker run -i mmsplice vep --plugin MMSplice --format vcf --assembly GRCh37 --database --port 3337 --vcf -o STDOUT | tee variant_effect_output.txt > $out_dir/working_files/$prefix.mmsplice.txt

    # Extract MMSplice scores
    ./extractMMSplice.sh -f $out_dir/working_files/$prefix.mmsplice.txt -p $prefix
    echo $(date +%x_%r) 'Finished running MMSplice'

    # Run SpliceAI
    echo $(date +%x_%r) 'Starting to run SpliceAI'

    spliceai -I $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf -D 1000 -R $reference_genome -A grch37 -O $out_dir/working_files/$prefix.spliceai.vcf
    ./extractSpliceAI.sh -f $out_dir/working_files/$prefix.spliceai.vcf
    echo $(date +%x_%r) 'Finished running SpliceAI'

fi

##################################
# STEP ?: Introme functions

echo $(date +%x_%r) 'Running AG Creation Check'
./AG_check.sh -p $prefix -r $reference_genome

echo $(date +%x_%r) 'Scoring insdels and MNVs'

## TO DO ##
## Implement a skip if no variants to score

echo $(bcftools view -h $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | grep -v "^#" | wc -l) 'MNV/insdel variants to score'
#./Testing/MNV.sh -r $reference_genome -p $prefix -f $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz

echo $(date +%x_%r) 'Second annotation'
vcfanno conf_introme.toml $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.scored.vcf.gz

##################################
# STEP 8: Output final list of variants as a spreadsheet-friendly TSV

if [ $no_csq == 0 ]; then
    bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/Variant_Type\t%INFO/Gene_Symbol\t%INFO/Gene_Strand\t%INFO/Gene_Location\t%INFO/Intron_Type\t%INFO/Gene_Regions\t%INFO/CADD_Phred\t%INFO/MGRB_AF\t%INFO/gnomAD_PM_AF\t%INFO/CSQ\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/MMSplice_Delta_Logit\t%INFO/AG_Created\t%INFO/GT_Created\t%INFO/Branchpointer_Branchpoint_Prob\t%INFO/Branchpointer_U2_Binding_Energy\t%INFO/SpliceAI_Distance\t%INFO/SpliceAI_Acceptor_Gain\t%INFO/SpliceAI_Acceptor_Loss\t%INFO/SpliceAI_Donor_Gain\t%INFO/SpliceAI_Donor_Loss\t%INFO/SpliceAI_DP_Acceptor_Gain\t%INFO/SpliceAI_DP_Acceptor_Loss\t%INFO/SpliceAI_DP_Donor_Gain\t%INFO/SpliceAI_DP_Donor_Loss[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.scored.vcf.gz > $out_dir/working_files/$prefix.introme.tsv
else
    bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/Variant_Type\t%INFO/Gene_Symbol\t%INFO/Gene_Strand\t%INFO/Gene_Location\t%INFO/Intron_Type\t%INFO/Gene_Regions\t%INFO/CADD_Phred\t%INFO/MGRB_AF\t%INFO/gnomAD_PM_AF\t.\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/MMSplice_Delta_Logit\t%INFO/AG_Created\t%INFO/GT_Created\t%INFO/Branchpointer_Branchpoint_Prob\t%INFO/Branchpointer_U2_Binding_Energy\t%INFO/SpliceAI_Distance\t%INFO/SpliceAI_Acceptor_Gain\t%INFO/SpliceAI_Acceptor_Loss\t%INFO/SpliceAI_Donor_Gain\t%INFO/SpliceAI_Donor_Loss\t%INFO/SpliceAI_DP_Acceptor_Gain\t%INFO/SpliceAI_DP_Acceptor_Loss\t%INFO/SpliceAI_DP_Donor_Gain\t%INFO/SpliceAI_DP_Donor_Loss[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.scored.vcf.gz > $out_dir/working_files/$prefix.introme.tsv
fi

# Remove the '[<number]' column numbers added by bcftools query to column names
grep '^#' $out_dir/working_files/$prefix.introme.tsv | sed "s/\[[0-9]*\]//g" > $out_dir/working_files/$prefix.introme.tsv.header
grep -v '^#' $out_dir/working_files/$prefix.introme.tsv >> $out_dir/working_files/$prefix.introme.tsv.header
mv $out_dir/working_files/$prefix.introme.tsv.header $out_dir/working_files/$prefix.annotated.tsv

# Sort by chromosome and coordinate
sort -k1,1n -k2,2n $out_dir/working_files/$prefix.annotated.tsv > $out_dir/working_files/$prefix.annotated.sorted.tsv
mv $out_dir/working_files/$prefix.annotated.sorted.tsv $out_dir/working_files/$prefix.annotated.tsv


##################################
# STEP ?: Generate ESE & ESS Scores
# python ESE_ESS_scoring.py $out_dir/working_files/$prefix.annotated.AG.tsv $out_dir/working_files/$prefix.annotated.AG.motifs.tsv $reference_genome

# STEP ?: Clean files
while read line; do
    if [[ ${line:0:1} == "#" ]]; then
        echo "$line" > $out_dir/$prefix.introme.tsv
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

        regions=$(echo "$line" | cut -f 11)
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

        echo "$line" | awk -v v="$variant_type" -v g="$gene_location" -v r="$region_replace" '{$7=v; $10=g; $11=r; print }' OFS='\t' >> $out_dir/$prefix.introme.tsv
    fi
done < $out_dir/working_files/$prefix.annotated.tsv

##################################
# STEP 6: Generate consensus scores

echo $(date +%x_%r) 'Beginning Introme consensus score generation'

# Run the appropriate machine learning model based off the method of score generation
if [ $mode == "full" ]; then
     Rscript --vanilla consensus_scoring.R full $out_dir/$prefix.introme.tsv $out_dir/$prefix.introme.predictions.tsv
elif [ $mode == "fast" ]; then
    #Remove commas from precomputed SpliceAI scores
    ./precomputed_processing.sh -f $out_dir/$prefix.introme.tsv -p $prefix
    mv $out_dir/$prefix.introme.nocommas.tsv $out_dir/$prefix.introme.tsv

     Rscript --vanilla consensus_scoring.R fast $out_dir/$prefix.introme.tsv $out_dir/$prefix.introme.fastpredictions.tsv
fi


##################################
# End

echo $(date +%x_%r) 'Introme complete'

if [ $mode == "full" ]; then
    echo "Your file can be found at $out_dir/$prefix.introme.predictions.tsv"
elif [ $mode == "fast" ]; then
    echo "Your file can be found at $out_dir/$prefix.introme.fastpredictions.tsv"
fi


exit
