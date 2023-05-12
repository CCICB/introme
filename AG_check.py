#!/usr/bin/env python3

import sys, shutil
import pysam

def main(input_file, reference_genome):
    # Open input VCF file, make copy of annotation VCF file
    input_file = pysam.VariantFile(input_file)
    shutil.copyfile("annotations/introme_annotate.vcf", "introme_annotate.functions2.vcf")

    # Open reference genome
    ref_genome = pysam.FastaFile(reference_genome)

    lines = []

    for i, record in enumerate(input_file):
        chrom = record.chrom
        pos = record.pos # VCF files are 1-based
        ref = record.ref
        alt = record.alts[0]
        strand = record.info.get("strand")

        if len(ref) == 1 and len(alt) == 1:
            variant_type = "SNV"
        elif len(ref) > 1 and len(alt) == 1:
            variant_type = "INDEL"
        elif len(ref) == 1 and len(alt) > 1:
            variant_type = "INDEL"
        elif len(ref) > 1 and len(alt) > 1:
            variant_type = "INSDEL"

        # start and end denote 0-based, half-open intervals.
        ref_seq = ref_genome.fetch(chrom, pos - 1 - 1, pos + len(ref))
        alt_seq = ref_seq[0] + alt + ref_seq[len(ref)+1:]

        # print(pos, ref_seq, alt_seq)

        pos_strand = strand.count("+") if strand is not None else 0
        neg_strand = strand.count("-") if strand is not None else 0

        append = ""

        if pos_strand > 0:
            if "AG" in alt_seq and "AG" not in ref_seq:
                append += "ag_created=+;"
            if "AG" not in alt_seq and "AG" in ref_seq:
                append += "ag_lost=+;"
            if "GT" in alt_seq and "GT" not in ref_seq:
                append += "gt_created=+;"
            if "GT" not in alt_seq and "GT" in ref_seq:
                append += "gt_lost=+;"

        if neg_strand > 0:
            ref_seq_rev = ref_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            alt_seq_rev = alt_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            if "AG" in alt_seq_rev and "AG" not in ref_seq_rev:
                append += "ag_created=-;"
            if "AG" not in alt_seq_rev and "AG" in ref_seq_rev:
                append += "ag_lost=-;"
            if "GT" in alt_seq_rev and "GT" not in ref_seq_rev:
                append += "gt_created=-;"
            if "GT" not in alt_seq_rev and "GT" in ref_seq_rev:
                append += "gt_lost=-;"

        lines.append("\t".join(str(item) for item in [chrom, pos, ".", ref, alt, ".", "PASS", f"variant_type={variant_type};{append}"]))

    with open('introme_annotate.functions2.vcf', 'a') as f:
        for line in lines:
            f.write(f"{line}\n")

if __name__ == "__main__":
    input_file = sys.argv[1]
    reference_genome = sys.argv[2]
    main(input_file, reference_genome)
