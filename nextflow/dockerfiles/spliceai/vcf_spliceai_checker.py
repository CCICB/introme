#!/usr/bin/env python3
import sys
import argparse

# List of field names corresponding to each transcript field.
FIELD_NAMES = [
    "ALLELE", "SYMBOL", "STRAND",
    "DS_AG", "DS_AL", "DS_DG", "DS_DL",
    "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    "DS_AG_REF", "DS_AG_ALT", "DS_AL_REF", "DS_AL_ALT",
    "DS_DG_REF", "DS_DG_ALT", "DS_DL_REF", "DS_DL_ALT"
]

def parse_spliceai(info_field):
    """
    Extract the SpliceAI value from the INFO field.
    """
    for entry in info_field.split(';'):
        if entry.startswith("SpliceAI="):
            return entry[len("SpliceAI="):]
    return None

def parse_transcript(transcript_str):
    """
    Parse a transcript string based on the expected 19 fields.
    Returns a dictionary mapping field names to their values.
    """
    fields = transcript_str.split('|')
    if len(fields) != len(FIELD_NAMES):
        return None, f"Expected {len(FIELD_NAMES)} fields but got {len(fields)} in transcript: {transcript_str}"

    transcript = {name: value for name, value in zip(FIELD_NAMES, fields)}
    return transcript, None

def compare_numeric(val1, val2):
    """
    Compare two numeric values provided as strings.
    If the value is a float (contains a decimal or exponent) compare within a tolerance of 0.01.
    Otherwise, compare as integers for exact equality.
    """
    try:
        if ('.' in val1 or 'e' in val1.lower()) or ('.' in val2 or 'e' in val2.lower()):
            f1, f2 = float(val1), float(val2)
            return abs(f1 - f2) <= 0.01 + 0.001
        else:
            return int(val1) == int(val2)
    except ValueError:
        return val1 == val2

def compare_transcript_entries(entry1, entry2):
    """
    Compare two transcript dictionaries field by field.
    For ALLELE and SYMBOL, values must match exactly.
    For all numeric fields, if they represent integers, they must be exactly equal; if floats, they must match within 0.01.
    Returns a list of differences with field names.
    """
    differences = []
    # Compare ALLELE and SYMBOL exactly
    for key in ["ALLELE", "SYMBOL"]:
        if entry1[key] != entry2[key]:
            differences.append(f"{key} differs: {entry1[key]} vs {entry2[key]}")

    # Compare all remaining fields numerically
    for key in FIELD_NAMES[2:]:
        if not compare_numeric(entry1[key], entry2[key]):
            differences.append(f"{key} differs: {entry1[key]} vs {entry2[key]}")
    return differences

def main():
    parser = argparse.ArgumentParser(
        description="Compare two VCF files for SpliceAI transcripts and scores (with field names) without exiting on first error."
    )
    parser.add_argument("vcf1", help="Path to the first VCF file")
    parser.add_argument("vcf2", help="Path to the second VCF file")
    args = parser.parse_args()

    errors = []

    # Read non-header variant lines from both VCF files.
    with open(args.vcf1, 'r') as f1, open(args.vcf2, 'r') as f2:
        lines1 = [line.strip() for line in f1 if not line.startswith('#') and line.strip()]
        lines2 = [line.strip() for line in f2 if not line.startswith('#') and line.strip()]

    if len(lines1) != len(lines2):
        errors.append(f"VCF files have different number of variants: {len(lines1)} vs {len(lines2)}")

    # Process up to the minimum number of variants.
    variant_count = min(len(lines1), len(lines2))
    for line_num in range(variant_count)[:100]:
        line1 = lines1[line_num]
        line2 = lines2[line_num]
        fields1 = line1.split('\t')
        fields2 = line2.split('\t')

        # Compare variant identity: chromosome, position, reference, and alternate.
        if (fields1[0] != fields2[0] or
            fields1[1] != fields2[1] or
            fields1[3] != fields2[3] or
            fields1[4] != fields2[4]):
            errors.append(f"Line {line_num+1}: Variant mismatch: {fields1[:5]} vs {fields2[:5]}")
            continue

        # Extract SpliceAI field from INFO column (8th column, index 7).
        spliceai1 = parse_spliceai(fields1[7])
        spliceai2 = parse_spliceai(fields2[7])
        if spliceai1 is None or spliceai2 is None:
            errors.append(f"Line {line_num+1}: SpliceAI key not found in INFO column")
            continue

        # Split the SpliceAI value into transcript entries.
        transcripts_str1 = spliceai1.split(',')
        transcripts_str2 = spliceai2.split(',')

        transcript_entries1 = []
        transcript_entries2 = []
        transcript_error = False

        for t_str in transcripts_str1:
            entry, err = parse_transcript(t_str)
            if err:
                errors.append(f"Line {line_num+1}: {err}")
                transcript_error = True
                break
            transcript_entries1.append(entry)
        for t_str in transcripts_str2:
            entry, err = parse_transcript(t_str)
            if err:
                errors.append(f"Line {line_num+1}: {err}")
                transcript_error = True
                break
            transcript_entries2.append(entry)
        if transcript_error:
            continue

        # Compare that both have the same number of transcript entries.
        if len(transcript_entries1) != len(transcript_entries2):
            errors.append(f"Line {line_num+1}: Number of transcript entries differ: {len(transcript_entries1)} vs {len(transcript_entries2)}")
            continue

        # Compare transcripts pairwise. We assume order is the same.
        for idx, (entry1, entry2) in enumerate(zip(transcript_entries1, transcript_entries2)):
            # For reporting, include SYMBOL to identify the transcript.
            identifier = entry1.get("SYMBOL", f"transcript {idx+1}")
            diffs = compare_transcript_entries(entry1, entry2)
            if diffs:
                errors.append(f"Line {line_num+1} {fields1[0]}/{fields1[1]}/{fields1[3]}/{fields1[4]} Transcript {identifier}: " + "; ".join(diffs))

    if errors:
        print("Differences found:")
        for err in errors:
            print(err)
    else:
        print("All variants and SpliceAI scores match.")

if __name__ == "__main__":
    main()
