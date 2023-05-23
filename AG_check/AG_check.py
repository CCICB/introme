#!/usr/bin/env python3

import subprocess, sys, shutil
import pysam
from Bio.Seq import Seq

def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) > 1 and len(alt) == 1:
        return "INDEL"
    elif len(ref) == 1 and len(alt) > 1:
        return "INDEL"
    elif len(ref) > 1 and len(alt) > 1:
        return "INSDEL"

def ag_gt_check(strand, alt_seq, ref_seq):
    pos_strand = 0 if strand is None else sum(value.count("+") for value in strand)
    neg_strand = 0 if strand is None else sum(value.count("-") for value in strand)

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
        ref_seq_rev = Seq(ref_seq).reverse_complement()
        alt_seq_rev = Seq(alt_seq).reverse_complement()
        
        if "AG" in alt_seq_rev and "AG" not in ref_seq_rev:
            append += "ag_created=-;"
        if "AG" not in alt_seq_rev and "AG" in ref_seq_rev:
            append += "ag_lost=-;"
        if "GT" in alt_seq_rev and "GT" not in ref_seq_rev:
            append += "gt_created=-;"
        if "GT" not in alt_seq_rev and "GT" in ref_seq_rev:
            append += "gt_lost=-;"

    return append

def main(input_vcf: pysam.VariantFile, reference_genome: pysam.FastaFile):
    # Make copy of introme annotation VCF header
    # TODO: use tmp files
    shutil.copyfile("annotations/introme_annotate.vcf", "introme_annotate.functions2.vcf")
    # vcf_out = pysam.VariantFile("introme_annotate.functions2.vcf", 'w', header=pysam.VariantFile("annotations/introme_annotate.vcf").header)

    with open('introme_annotate.functions2.vcf', 'a') as f:
        for record in input_vcf:
            chrom = record.chrom
            pos = record.pos # VCF files are 1-based
            ref = record.ref
            alt = record.alts[0]
            strand = record.info.get("strand")

            variant = variant_type(ref, alt)

            # In pysam.FastaFile.fetch(), start and end denote 0-based, half-open intervals.
            ref_seq = reference_genome.fetch(chrom, pos - 1 - 1, pos + len(ref))
            alt_seq = ref_seq[0] + alt + ref_seq[len(ref) + 1:]

            # print(chrom, pos, ref, alt, strand)
            # print(pos, ref_seq, alt_seq)

            append = ag_gt_check(strand, alt_seq, ref_seq)

            f.write("\t".join(str(item) for item in [chrom, pos, ".", ref, alt, ".", "PASS", f"variant_type={variant};{append}"]) + "\n")
            
            # vcf_out.new_record(chrom=chrom, pos=pos, ref=ref, alt=alt, filter="PASS", info=f"variant_type={variant_type};{append}")
            # f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tvariant_type={variant};{append}\n")
            # ("\t".join(str(item) for item in [chrom, pos, ".", ref, alt, ".", "PASS", f"variant_type={variant_type};{append}"]))

    # # bcftools sort introme_annotate.functions.vcf
    # sort_process = subprocess.Popen(['bcftools', 'sort', 'introme_annotate.functions2.vcf'], stdout=subprocess.PIPE)

    # # uniq
    # uniq_process = subprocess.Popen(['uniq'], stdin=sort_process.stdout, stdout=subprocess.PIPE)
    # sort_process.stdout.close()

    # # bgzip > introme_annotate.functions.vcf.gz
    # bgzip_process = subprocess.Popen(['bgzip'], stdin=uniq_process.stdout, stdout=open('introme_annotate.functions2.vcf.gz', 'wb'))
    # uniq_process.stdout.close()
    # bgzip_process.wait()

    # # rm introme_annotate.functions.vcf
    # subprocess.run(['rm', 'introme_annotate.functions2.vcf'])

    # # tabix -f introme_annotate.functions.vcf.gz
    # subprocess.run(['tabix', '-f', 'introme_annotate.functions2.vcf.gz'])


if __name__ == "__main__":
    input_vcf = pysam.VariantFile(sys.argv[1])
    reference_genome = pysam.FastaFile(sys.argv[2])

    main(input_vcf, reference_genome)
