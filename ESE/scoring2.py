from motifs import RBPsplice, RBPmotifs
from variants import Variant, VariantContext

import pysam
import csv
import sys
import numpy as np

CONTEXT_LENGTH = 7

def read_vcf(reader: csv.reader, writer: csv.writer, reference_genome: pysam.FastaFile, RBPmotifs: list[RBPsplice]):
    for row in reader:
        chromosome = str(row[0])

        # Write header line
        if chromosome[:2] == "##":
            continue
        elif chromosome == "#CHROM":
            row.extend([motif.name for motif in RBPmotifs])
            writer.writerow(row)
            continue

        position = int(row[1])
        ref = str(row[3])
        alt = str(row[4])
        strand = "+" #TODO: get strand; deal with double strand variants

        # Extract reference sequence
        # Skip insertions longer than 30bp
        if len(alt) > 30:
            row.extend(['.','.','.','.','.','.'])
            writer.writerow(row)
            continue
        
        variant = Variant("chr" + chromosome, position, ref, alt)
        variant_context = variant.faidx_context(reference_genome, CONTEXT_LENGTH)

        for motif in RBPmotifs:
            # print(motif.name)
            # print(motif.calculate(variant_context.ref_sequence(motif.length)))
            # print(motif.calculate(variant_context.alt_sequence(motif.length)))

             a = round(
                np.max(motif.calculate(variant_context.alt_sequence(motif.length)))
                - np.max(motif.calculate(variant_context.ref_sequence(motif.length)))
            , 2)
             
             a = "0" if a == 0 else a
             row.append(a)

        writer.writerow(row)

def main():
    CONTEXT_LENGTH = 7
    for motif in RBPmotifs:
        print(motif.name, motif.threshold)
        print(motif)

    file = sys.argv[1]
    output = sys.argv[2]
    reference_genome = pysam.FastaFile(sys.argv[3])

    with (
        pysam.FastaFile(sys.argv[3]) as reference_genome,
        open(file) as read, open(output, 'w+') as write
    ):
        reader = csv.reader(read, delimiter='\t')
        writer = csv.writer(write, delimiter='\t')

        read_vcf(reader, writer, reference_genome, RBPmotifs)


if __name__ == "__main__":
    main()
