# /bin/sh
# In case you download dbscSNV and it comes separately per chromosome

# Initialize the combined file with the header from chr1
head -n 1 dbscSNV1.1.chr1 > combined_annotations.tsv

# Iterate through chromosomes 1-22, then X and Y
for chr in {1..22} X Y; do
    # Check if the file exists to avoid errors
    if [[ -f "dbscSNV1.1.chr${chr}" ]]; then
        # Append data excluding the header
        tail -n +2 "dbscSNV1.1.chr${chr}" >> combined_annotations.tsv
    else
        echo "Warning: dbscSNV1.1.chr${chr} does not exist and will be skipped."
    fi
done

echo "Merging complete. Combined file: combined_annotations_sorted.tsv"

