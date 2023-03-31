#!/bin/bash

# Introme Bash Script
#
# Developers:
# Patricia Sullivan (psullivan@ccia.org.au)
# Velimir Gayevskiy (vel@vel.nz)
# Sarah Beecroft (sarah.beecroft@uwa.edu.au)


# Version=1.0.0
##################################
# DEFINE VARIABLES

# Directories
out_dir='./output'

# Trigger values
help=0
no_pop_filter=0
no_qual_filter=0
toggle_safety=1
mode=full

# Hard filter cutoffs for each annotation, filtering thresholds will depend on your application, we urge you to carefully consider these
max_AF='0.01' # Minimum gnomAD allele frequency in the population with the highest allele frequency in gnomAD
min_QUAL='>=200' # The QUAL VCF field
min_DP='>=20' # The sample with the highest depth must have a equal or larger depth than this value
min_AD='>=5' # The sample with the highest number of alternate reads must have a equal or larger number than this value

##################################
# INTROME ARGUMENTS

while getopts "a:b:f:g:hm:p:qr:sv:" opt; do
    case $opt in
        a) genome="$OPTARG";; # Genome assembly (hg19 or hg38)
        b) input_BED="$OPTARG";; # Input BED file (i.e. regions of interest)
        f) max_AF="$OPTARG";; # maximum allele frequency allowed for variants
        g) input_gtf="$OPTARG";; # Path to the gtf used for MMSplice
        h) help=1;; # trigger help command
        m) mode="$OPTARG";; # Mode in which introme is running - add "fast" for precomputed scores only
        p) prefix="$OPTARG";; # Output file prefix
        q) no_qual_filter=1;; # trigger no quality filter command
        r) reference_genome="$OPTARG";; # Path to the reference genome used for mapping
        s) toggle_safety=0;; # turn off Introme single score check
        v) input_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

# Run help command
if [[ $help == 1 ]]; then
    echo "Program: Introme (A program for variant calling and manipulating VCFs and BCFs)"
    echo "Version: 1.0.0"$'\n'

    echo "Usage: ./run_introme.sh [options] -b <subset.bed.gz> -p <prefix> -r <reference_genome.fa> -v <variants.vcf.gz>"$'\n'

    echo "Required Commands:"
    echo "-g <annotations.gtf>        Input GTF file (protein-coding only preferred)"
    echo "-p <prefix>                 Output file prefix"
    echo "-r <reference_genome.fa>    The reference genome used for mapping"
    echo "-v <variants.vcf.gz>        Input variants in VCF format (must be gzipped)"$'\n'

    echo "Options:"
    echo "-a <hg19/hg38>              Genome assembly used in VCF"
    echo "-b <subset.bed.gz>          Input BED file of the regions of interest"
    echo "-f <allele_frequency>       Maximum allele frequency to include"
    echo "-h                          Print Introme usage instructions"
    echo "-q                          Remove variant quality filtering"
    echo "-s                          Turn off Introme single score check"$'\n'
    exit 1
fi

# Make sure each required argument has a non-empty value
if [ -z $input_gtf ]; then
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
elif [ -z $genome ]; then
    if [[ $(echo "$reference_genome" | grep -c -e "38") > 0 ]]; then
        genome='hg38'
    elif [[ $(echo "$reference_genome" | grep -c -e "19" -e "37") > 0 ]]; then
        genome='hg19'
    else
        echo "Genome assembly could not be determined. Please add either '-a hg19' or '-a hg38'"
    	exit 1
    fi
fi

echo "$genome"


##################################
# STEP 1: subsetting the VCF to genomic regions of interest (first because it gets rid of the most variants)

echo $(date +%x_%r) 'Beginning subsetting to genomic regions of interest'
echo $(gzip -d -c $input_VCF | grep -v '^#' | wc -l) 'variants prior to subsetting'

bcftools sort $input_VCF | uniq | bgzip > $out_dir/working_files/$prefix.sorted.vcf.gz  # Ensures the file is sorted correctly prior to subsetting
bcftools norm -m-both $out_dir/working_files/$prefix.sorted.vcf.gz | bgzip > $out_dir/working_files/$prefix.sorted.norm.vcf.gz

if [ -z $input_BED ]; then
    echo $(date +%x_%r) 'No BED file provided - Beginning subsetting to GTF regions'
    bedtools intersect -header -u -a $out_dir/working_files/$prefix.sorted.norm.vcf.gz -b $input_gtf | bgzip > $out_dir/working_files/$prefix.subset.vcf.gz # -u for unique record in VCF
else
    echo $(date +%x_%r) 'BED file provided - Beginning subsetting to genomic regions of interest'
    bedtools intersect -header -u -a $out_dir/working_files/$prefix.sorted.norm.vcf.gz -b $input_BED | bgzip > $out_dir/working_files/$prefix.subset.vcf.gz # -u for unique record in VCF
fi

tabix -p vcf $out_dir/working_files/$prefix.subset.vcf.gz

variant_count=$(bcftools view -H $out_dir/working_files/$prefix.subset.vcf.gz | wc -l | tr -d ' ')
echo $(date +%x_%r) 'Subsetting complete -' $variant_count 'variants remaining'

if [[ $variant_count == 0 ]]; then
    exit 1
fi

##################################
# STEP 2: Hard filtering on variant quality (this is here to reduce the number of variants going into the CPU-costly annotation step below)

echo $(date +%x_%r) 'Beginning quality filtering'

if [[ $no_qual_filter == 0 ]]; then
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.') && MAX(FORMAT/DP[*])$min_DP && MAX(FORMAT/AD[*:1])$min_AD" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
else
    bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(FILTER='PASS' || FILTER='.') && (QUAL$min_QUAL || QUAL='.')" $out_dir/working_files/$prefix.subset.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.vcf.gz
fi

tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.vcf.gz

variant_count=$(bcftools view -H $out_dir/working_files/$prefix.subset.highquality.vcf.gz | wc -l | tr -d ' ')
echo $(date +%x_%r) 'Quality filtering complete -' $variant_count 'variants remaining'

if [[ $variant_count == 0 ]]; then
    exit 1
fi

##################################
# STEP 3: annotate the subsetted VCF with useful information, to be used for filtering downstream


export IRELATE_MAX_GAP=1000 # Note: this is set to speed up annotation when .csi (as opposed to .tbi) files are present via https://github.com/brentp/vcfanno/releases/

# Run the appropriate annotation command based off the mode supplied
if [ $mode == "full" ]; then
    echo $(date +%x_%r) 'Beginning annotation'
    vcfanno -p $(getconf _NPROCESSORS_ONLN) -lua conf.lua annotations/gencode.$genome.toml $out_dir/working_files/$prefix.subset.highquality.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz # The conf.toml file specifies what VCFanno should annotate
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
# STEP 4: Hard filtering on the values of annotations added in the previous step

echo $(date +%x_%r) 'Beginning filtering by annotation values'

bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i"(gnomAD_PM_AF<=$max_AF || gnomAD_PM_AF='.')" $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz
tabix -p vcf $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz

variant_count=$(bcftools view -H $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | wc -l | tr -d ' ')
echo $(date +%x_%r) 'Frequency filtering complete -' $variant_count 'variants remaining'

if [[ $variant_count == 0 ]]; then
    exit 1
fi

##################################
# STEP 5: Run MMSplice and Splice AI - Skipped if fast mode

if [ $mode == "full" ]; then
    gunzip -k $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz # MMSplice and SpliceAI needs unzipped input files

    # Remove annotations which interfere with score extraction
    bcftools view -h $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz > $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
    grep -v "^#" $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf | awk '{$8="."; print }' OFS='\t' >> $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
    bgzip -f $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf
    gunzip -k $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf

    # Run MMSplice using docker
    echo $(date +%x_%r) 'Starting to run MMSplice'
    docker run -v $(pwd):/data mmsplice \
        --vcf /data/$out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf.gz \
        --gtf /data/$input_gtf --output /data/$out_dir/working_files/$prefix.mmsplice.vcf \
        --fasta /data/$reference_genome

    # Extract MMSplice scores
    cp annotations/introme_annotate.vcf introme_annotate.mmsplice.vcf
    # Continues process if no variants were scored by MMSplice
    if [[ -z $(ls $out_dir/working_files/$prefix.mmsplice.vcf) ]]; then
        cp $out_dir/working_files/test.mmsplice.vcf $out_dir/working_files/$prefix.mmsplice.vcf
    fi

    bcftools view -H $out_dir/working_files/$prefix.mmsplice.vcf | cut -f1-8 | grep "CSQ" | \
        sed -e 's/CSQ=/alt_acceptor=/' -e 's/|/;alt_acceptor_intron=/' -e 's/|/;alt_donor=/' \
        -e 's/|/;alt_donor_intron=/' -e 's/|/;alt_exon=/' -e 's/|/;delta_logit_PSI=/' \
        -e 's/|/;pathogenicity=/' -e 's/|/;ref_acceptor=/' -e 's/|/;ref_acceptor_intron=/' \
        -e 's/|/;ref_donor=/' -e 's/|/;ref_donor_intron=/' -e 's/|/;ref_exon=/' \
        >> introme_annotate.mmsplice.vcf
    bcftools sort introme_annotate.mmsplice.vcf -O z -o introme_annotate.mmsplice.vcf.gz
    tabix introme_annotate.mmsplice.vcf.gz

    echo $(date +%x_%r) 'Finished running MMSplice'

    # # Run SpliceAI
    echo $(date +%x_%r) 'Starting to run SpliceAI'

    spliceai -I $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf -D 1000 -R $reference_genome -A grch37 -O $out_dir/working_files/$prefix.spliceai.vcf
    ./extractSpliceAI.sh -p $prefix

    echo $(date +%x_%r) 'Finished running SpliceAI'

    echo $(date +%x_%r) 'Starting to run Spliceogen'
    docker run -v $(pwd):/data spliceogen \
        -input /data/$out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf \
        -fasta /data/$reference_genome \
        -gtf /data/$input_gtf

    # Extract Spliceogen Scores
    cp annotations/introme_annotate.vcf introme_annotate.spliceogen.vcf
    grep -v "^#" $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf_out.txt | \
        cut -f1,2,4,5,8-11,16-23 | sed -e $'s/\t/\t.\t/2' -e $'s/\t/\t.\tPASS\t/5' \
        -e $'s/\t/\tmesDonRef=/7' -e $'s/\t/;mesDonAlt=/8' -e $'s/\t/;mesAccRef=/8' \
        -e $'s/\t/;mesAccAlt=/8' -e $'s/\t/;ESEmaxRef=/8' -e $'s/\t/;ESEmaxAlt=/8' \
        -e $'s/\t/;ESEminRef=/8' -e $'s/\t/;ESEminAlt=/8' -e $'s/\t/;DonGainP=/8' \
        -e $'s/\t/;AccGainP=/8' -e $'s/\t/;DonLossP=/8' -e $'s/\t/;AccLossP=/8' \
        >> introme_annotate.spliceogen.vcf
    bcftools sort introme_annotate.spliceogen.vcf -O z -o introme_annotate.spliceogen.vcf.gz
    tabix introme_annotate.spliceogen.vcf.gz

    echo $(date +%x_%r) 'Finished running Spliceogen'

fi

#################################
# STEP 6: Introme functions

echo $(date +%x_%r) 'Running AG Creation Check'
./AG_check.sh -i $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz -r $reference_genome

echo $(date +%x_%r) 'Running ESE Scoring'
python ESE_ESS_scoring.py $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.rmanno.vcf $out_dir/working_files/$prefix.subset.highquality.ESE.vcf $reference_genome
cat annotations/introme_annotate.vcf $out_dir/working_files/$prefix.subset.highquality.ESE.vcf | bgzip > introme_annotate.ESE.vcf

bcftools sort introme_annotate.ESE.vcf | uniq | bgzip > introme_annotate.ESE.vcf.gz
rm introme_annotate.ESE.vcf
tabix -f introme_annotate.ESE.vcf.gz

echo $(date +%x_%r) 'Scoring insdels and MNVs'

MNVs=$(bcftools filter -i"TYPE!='snp' && TYPE!='indel'" $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz | grep -v "^#" | wc -l | tr -d ' ')

if [[ $MNVs > 0 ]]; then
    echo $MNVs 'MNV/insdel variants to score'
    ./Testing/MNV.sh -r $reference_genome -p $prefix -f $out_dir/working_files/$prefix.subset.highquality.annotated.filtered.vcf.gz
else
    echo 'No MNV/insdel variants to score'
fi

echo $(date +%x_%r) 'Splicing annotation'
vcfanno -p $(getconf _NPROCESSORS_ONLN) -lua conf.lua annotations/annotate.$genome.toml $out_dir/working_files/$prefix.subset.highquality.annotated.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.splicing_anno.vcf.gz
tabix -f $out_dir/working_files/$prefix.subset.highquality.annotated.splicing_anno.vcf.gz

vcfanno -p $(getconf _NPROCESSORS_ONLN) -lua conf.lua annotations/conf_introme.toml $out_dir/working_files/$prefix.subset.highquality.annotated.splicing_anno.vcf.gz | bgzip > $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz
tabix -f $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz

##################################
# STEP 7: Output final list of variants as a spreadsheet-friendly TSV

# Check if CSQ column is present (VEP)
if [[ ! -z $(bcftools view $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz | grep "CSQ") ]]; then
    bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/Variant_Type\t%INFO/gene\t%INFO/strand\t%INFO/Gene_Location\t%INFO/Intron_Type\t%INFO/U12_score\t%INFO/Gene_Regions\t%INFO/CSQ\t%INFO/CADD_Phred\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/mesAccRef\t%INFO/mesAccAlt\t%INFO/mesDonRef\t%INFO/mesDonAlt\t%INFO/AG_Created\t%INFO/AG_Lost\t%INFO/GT_Lost\t%INFO/GT_Created\t%INFO/Branchpointer_Prob\t%INFO/Branchpointer_U2_Binding_Energy\t%INFO/SpliceAI_Acceptor_Gain\t%INFO/SpliceAI_Acceptor_Loss\t%INFO/SpliceAI_Donor_Gain\t%INFO/SpliceAI_Donor_Loss\t%INFO/DP_AG\t%INFO/DP_AL\t%INFO/DP_DG\t%INFO/DP_DL\t%INFO/AccGainP\t%INFO/AccLossP\t%INFO/DonGainP\t%INFO/DonLossP\t%INFO/MMSplice_alt_acceptor\t%INFO/MMSplice_alt_acceptor_intron\t%INFO/MMSplice_alt_donor\t%INFO/MMSplice_alt_donor_intron\t%INFO/MMSplice_alt_exon\t%INFO/MMSplice_delta_logit_PSI\t%INFO/MMSplice_pathogenicity\t%INFO/MMSplice_ref_acceptor\t%INFO/MMSplice_ref_acceptor_intron\t%INFO/MMSplice_ref_donor\t%INFO/MMSplice_ref_donor_intron\t%INFO/MMSplice_ref_exon\t%INFO/SRSF1_ref\t%INFO/SRSF1_alt\t%INFO/SRSF1_igM_ref\t%INFO/SRSF1_igM_alt\t%INFO/SRSF2_ref\t%INFO/SRSF2_alt\t%INFO/SRSF5_ref\t%INFO/SRSF5_alt\t%INFO/SRSF6_ref\t%INFO/SRSF6_alt\t%INFO/hnRNPA1_ref\t%INFO/hnRNPA1_alt\t%INFO/ESEmaxRef\t%INFO/ESEmaxAlt\t%INFO/ESEminRef\t%INFO/ESEminAlt[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz > $out_dir/working_files/$prefix.annotated.tsv
else
    bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/Variant_Type\t%INFO/gene\t%INFO/strand\t%INFO/Gene_Location\t%INFO/Intron_Type\t%INFO/U12_score\t%INFO/Gene_Regions\t%INFO/CADD_Phred\t%INFO/SPIDEX_dPSI_Max_Tissue\t%INFO/SPIDEX_dPSI_Zscore\t%INFO/dbscSNV_AdaBoost_Score\t%INFO/dbscSNV_RandomForest_Score\t%INFO/mesAccRef\t%INFO/mesAccAlt\t%INFO/mesDonRef\t%INFO/mesDonAlt\t%INFO/AG_Created\t%INFO/AG_Lost\t%INFO/GT_Lost\t%INFO/GT_Created\t%INFO/Branchpointer_Prob\t%INFO/Branchpointer_U2_Binding_Energy\t%INFO/SpliceAI_Acceptor_Gain\t%INFO/SpliceAI_Acceptor_Loss\t%INFO/SpliceAI_Donor_Gain\t%INFO/SpliceAI_Donor_Loss\t%INFO/DP_AG\t%INFO/DP_AL\t%INFO/DP_DG\t%INFO/DP_DL\t%INFO/AccGainP\t%INFO/AccLossP\t%INFO/DonGainP\t%INFO/DonLossP\t%INFO/MMSplice_alt_acceptor\t%INFO/MMSplice_alt_acceptor_intron\t%INFO/MMSplice_alt_donor\t%INFO/MMSplice_alt_donor_intron\t%INFO/MMSplice_alt_exon\t%INFO/MMSplice_delta_logit_PSI\t%INFO/MMSplice_pathogenicity\t%INFO/MMSplice_ref_acceptor\t%INFO/MMSplice_ref_acceptor_intron\t%INFO/MMSplice_ref_donor\t%INFO/MMSplice_ref_donor_intron\t%INFO/MMSplice_ref_exon\t%INFO/SRSF1_ref\t%INFO/SRSF1_alt\t%INFO/SRSF1_igM_ref\t%INFO/SRSF1_igM_alt\t%INFO/SRSF2_ref\t%INFO/SRSF2_alt\t%INFO/SRSF5_ref\t%INFO/SRSF5_alt\t%INFO/SRSF6_ref\t%INFO/SRSF6_alt\t%INFO/hnRNPA1_ref\t%INFO/hnRNPA1_alt\t%INFO/ESEmaxRef\t%INFO/ESEmaxAlt\t%INFO/ESEminRef\t%INFO/ESEminAlt[\t%GT][\t%DP][\t%AD][\t%GQ]\n' $out_dir/working_files/$prefix.subset.highquality.annotated.scored.vcf.gz > $out_dir/working_files/$prefix.annotated.tsv
fi

# Remove the '[<number]' column numbers added by bcftools query to column names
grep '^#' $out_dir/working_files/$prefix.annotated.tsv | sed "s/\[[0-9]*\]//g" > $out_dir/working_files/$prefix.annotated.format.tsv
grep -v '^#' $out_dir/working_files/$prefix.annotated.tsv >> $out_dir/working_files/$prefix.annotated.format.tsv

# Sort by chromosome and coordinate
sort -k1,1n -k2,2n $out_dir/working_files/$prefix.annotated.format.tsv > $out_dir/working_files/$prefix.annotated.tsv
rm $out_dir/working_files/$prefix.annotated.format.tsv

##################################
# STEP 8: Generate consensus scores

echo $(date +%x_%r) 'Beginning Introme consensus score generation'

#Run the appropriate machine learning model based off the method of score generation
if [ $mode == "full" ]; then
     Rscript --vanilla consensus_scoring.R full $out_dir/working_files/$prefix.annotated.tsv $out_dir/$prefix.introme.tsv
     echo "Your file can be found at $out_dir/$prefix.introme.predictions.tsv"
elif [ $mode == "fast" ]; then
    Rscript --vanilla consensus_scoring.R fast $out_dir/working_files/$prefix.annotated.tsv $out_dir/$prefix.introme.fastpredictions.tsv
    echo "Your file can be found at $out_dir/$prefix.introme.fastpredictions.tsv"
fi

##################################
# End

echo $(date +%x_%r) 'Introme complete'

exit
