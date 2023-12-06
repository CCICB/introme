# A file to update the input values of the params.json file
# So that the introme nextflow pipeline can run on user provided files

# Introme Arguments
while getopts "a:b:f:g:h:p:qr:sv:" opt; do
    case $opt in
        a) genome="$OPTARG";; # Genome assembly (hg19 or hg38)
        b) input_BED="$OPTARG";; # Input BED file (i.e. regions of interest)
        f) max_AF="$OPTARG";; # maximum allele frequency allowed for variants
        g) input_gtf="$OPTARG";; # Path to the gtf used for MMSplice
        h) help=1;; # trigger help command
        p) prefix="$OPTARG";; # Output file prefix
        q) no_qual_filter=1;; # trigger no quality filter command
        r) reference_genome="$OPTARG";; # Path to the reference genome used for mapping
        v) input_VCF="$OPTARG";; # Input VCF file
    esac
done
shift $(( $OPTIND - 1 )) # Remove parsed options and args from $@ list

# Run help command
if [[ $help == 1 ]]; then
    echo "Usage: ./update_params.sh [options] -b <subset.bed.gz> -p <prefix> -r <reference_genome.fa> -v <variants.vcf.gz> -g <annotations.gtf>"$'\n'

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

# If not in the direcoty download the params file - otherwise unnecessary
# wget https://raw.githubusercontent.com/CCICB/introme/nextflow/nextflow/params.json

# Make temporary version of params file to make edits to
cp params.json temp_params.json

# Make changes to params file based on the arguments and options provided
# Required Arguments
if [ -n $input_gtf ]; then
    sed "s|input/gtf|$input_gtf|" params.json > temp_params.json
    mv temp_params.json params.json
fi

if [ -n $input_VCF ]; then
    sed "s|input/toy.vcf.gz|$input_VCF|" params.json > temp_params.json
    mv temp_params.json params.json
fi

if [ -n $reference_genome ]; then
    sed "s|input/ref_genome|$reference_genome|" params.json > temp_params.json
    mv temp_params.json params.json
fi

if [ -n $genome ]; then
    sed "s|hg38|$genome|" params.json > temp_params.json
    mv temp_params.json params.json
fi

# Options
if [ ! -z $input_BED ]; then
    sed "s|input/toy.bed.gz|$input_BED|" params.json > temp_params.json
    mv temp_params.json params.json
fi

if [ ! -z $max_AF ]; then
    sed "s|0.01|$max_AF|" params.json > temp_params.json
    mv temp_params.json params.json
fi

if [ ! -z $prefix ]; then
    sed "s|splice_test_04_10|$prefix|" params.json > temp_params.json
    mv temp_params.json params.json
fi

if [ ! -z $no_qual_filter ]; then
    sed "s|true|$no_qual_filter|" params.json > temp_params.json
    mv temp_params.json params.json
fi

