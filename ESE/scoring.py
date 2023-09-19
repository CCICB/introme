from motifs import RBPsplice
from variants import Variant, StrandDirection
from ESEfinder_motif_source import ESEfinder_motifs

import pysam
import csv
import sys
import os
from enum import Enum, auto
import numpy as np
import pandas as pd

CONTEXT_LENGTH = 14

strand_d = {
    'strand=+': StrandDirection.FORWARD,
    'strand=-': StrandDirection.REVERSE,
    'strand=+,-': StrandDirection.BOTH,
    '.': StrandDirection.UNKNOWN
}

def read_vcf(vcf: pysam.VariantFile, ref_genome: pysam.FastaFile, RBPmotifs: list[RBPsplice]) -> pd.DataFrame:
    data: list[list] = []

    for record in vcf:
        chromosome = record.chrom
        position = record.pos
        id_ = record.id
        ref = record.ref
        # alt
        quality = record.qual
        filter_ = record.filter
        info = record.info

        assert(record.alts is not None and len(record.alts) == 1), (
            f"{chromosome}:{position}-{ref}>{record.alts}. Records should be left normalised"
        )

        alt = record.alts[0]
        strand = ';'.join(str(s) for s in record.info['strand'])

        if not "+" in strand and not "-" in strand:
            strand_dir = StrandDirection.UNKNOWN
        elif "+" in strand:
            strand_dir = StrandDirection.FORWARD
        elif "-" in strand:
            strand_dir = StrandDirection.REVERSE
        else:
            strand_dir = StrandDirection.BOTH


        # TODO pandas df or something
        # # Write header line
        # if chromosome[:2] == "##":
        #     continue
        # elif chromosome == "#CHROM":
        #     row.extend([motif.name for motif in RBPmotifs])
        #     writer.writerow(row)
        #     continue

        # Extract reference sequence
        # Skip insertions longer than 50bp
        # if len(alt) > 50:
        #     row.extend(['.','.','.','.','.','.'])
        #     writer.writerow(row)
        #     continue
        
        variant = Variant("chr" + chromosome, position, "." if ref is None else ref, alt, strand_dir)
        variant_context = variant.faidx_context(ref_genome, CONTEXT_LENGTH)

        motif_scores = []
        for motif in RBPmotifs:
            # print(motif.name)
            # print(motif.calculate(variant_context.ref_sequence(motif.length)))
            # print(motif.calculate(variant_context.alt_sequence(motif.length)))

            a = motif.calculate_variant(variant_context)
            a = "0" if a == 0 else a
            motif_scores.append(a)
            # append motif score to row
            # row.append(a)
        
        df_row = [chromosome, position, id_, ref, alt, quality, filter_, info] + motif_scores
        data.append(df_row)

    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] + [motif.name for motif in RBPmotifs]
    df = pd.DataFrame(data, columns=columns)
    return df

def is_path_writable(path: str) -> bool:
    """Check if a file path is writable."""
    
    # Check if directory exists
    if not os.path.isdir(os.path.dirname(path)):
        return False
    
    # If file exists, check if it's writable
    if os.path.exists(path):
        return os.access(path, os.W_OK)
    
    # If file doesn't exist, try to create it to check writability
    try:
        open(path, 'a').close()   # open in append mode and immediately close
        os.remove(path)           # remove the file after the test
        return True
    except Exception:
        return False

def main():
    CONTEXT_LENGTH = 14
    motifs = []
    for _, motif in ESEfinder_motifs.motifs.items():
        motifs.append(RBPsplice.from_2D_list(motif['matrix'], motif['name'], threshold=motif['threshold']))

    vcf_file =  pysam.VariantFile(sys.argv[1])
    output_path = sys.argv[2]
    reference_genome = pysam.FastaFile(sys.argv[3])

    if not is_path_writable(output_path):
        raise ValueError(f"File path '{output_path}' is not writable!")

    df = read_vcf(vcf_file, reference_genome, motifs)

    df.to_csv(output_path, encoding='utf-8', index=False)

if __name__ == "__main__":
    main()
