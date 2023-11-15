import csv
import sys

def read_tsv(filename):
    variants = {}
    with open(filename, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fail = 0
        for i, row in enumerate(reader):
            key = ("chr" + row['chr'], row['pos_hg38'], row['ref'], row['alt'])
            if key in variants:
                print(f"row {i}: {key} is already seen")
                if variants[key] != row:
                    print(variants[key], row)
                    # exit(0)
                    fail += 1

            if fail > 20: break
            variants[key] = row

    return variants

def read_vcf(filename):
    variants = {}
    with open(filename, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', fieldnames=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
        fail = 0
        for i, row in enumerate(reader):
            if row['#CHROM'].startswith('#') or row['#CHROM'] == "CHROM":
                continue  # skip header lines
            info_fields = row['INFO'].split(';')

            strand =''
            for field in info_fields:
                if field.startswith('strand='):
                    strand = field.split('=')[1]

            assert(strand != '')

            key = (row['#CHROM'], row['POS'], row['REF'], row['ALT'])
            if key in variants:
                print(f"row {i}: {key} is already seen")
                if variants[key] != row:
                    print(variants[key], row)
                    # exit(0)
                    fail += 1

            if fail > 20: break

            variants[key] = row

    return variants

def compare_variants(tsv_variants, vcf_variants):
    overlapping_variants = {}
    unique_tsv_variants = {}
    unique_vcf_variants = {}
    
    for key in tsv_variants:
        if key in vcf_variants:
            overlapping_variants[key] = (tsv_variants[key], vcf_variants[key])
        else:
            unique_tsv_variants[key] = tsv_variants[key]
    
    for key in vcf_variants:
        if key not in tsv_variants:
            unique_vcf_variants[key] = vcf_variants[key]
    
    return overlapping_variants, unique_tsv_variants, unique_vcf_variants

# Paths to your files
tsv_file = sys.argv[1]
vcf_file = sys.argv[2]

tsv_variants = read_tsv(tsv_file)
vcf_variants = read_vcf(vcf_file)

overlapping_variants, unique_tsv_variants, unique_vcf_variants = compare_variants(tsv_variants, vcf_variants)

# Now you have three dictionaries:
# overlapping_variants contains the variants present in both files
# unique_tsv_variants contains the variants only present in the TSV file
# unique_vcf_variants contains the variants only present in the VCF file

print(len(overlapping_variants), len(unique_tsv_variants), len(unique_vcf_variants))

