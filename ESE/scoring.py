from motifs import RBPsplice
from variants import Variant, VariantContext, StrandDirection
from ESEfinder_motif_source import ESEfinder_motifs
from dataclasses import dataclass
from RCRUNCH_motif_source import RCRUNCH_motifs
from MSE import doMSE, diff_in_scoring
from typing import Iterator, TextIO, Optional

import pysam
import csv
import sys
import os
from enum import Enum, auto
import numpy as np
import pandas as pd

CONTEXT_LENGTH = 22

# strand_d = {
#     'strand=+': StrandDirection.FORWARD,
#     'strand=-': StrandDirection.REVERSE,
#     'strand=+,-': StrandDirection.BOTH,
#     '.': StrandDirection.UNKNOWN
# }

@dataclass(frozen=True)
class VcfInfo():
    chromosome: any
    position: any
    id_: any
    ref: any
    alt: any
    quality: any
    filter_: any
    info: any

    def toList(self) -> list:
        return [self.chromosome, self.position, self.id_, self.ref, self.alt, self.quality, self.filter_, self.info]

def calculate_variants(variants: Iterator[tuple[Variant, VcfInfo]], ref_genome: pysam.FastaFile, RBPmotifs: list[RBPsplice]) -> pd.DataFrame:
    data: list[list] = []

    for variant, vcf_info in variants:      
        variant_context = variant.faidx_context(ref_genome, CONTEXT_LENGTH)
        if variant_context is None:
            continue
        motif_scores = calculuate_motifs(RBPmotifs, variant_context)

        if motif_scores is None:
            continue
        
        df_row = vcf_info.toList() + motif_scores
        data.append(df_row)

    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    expanded_motif_names = []
    for motif_name in [motif.name for motif in RBPmotifs]:
        expanded_motif_names.append(f"{motif_name}_ref")
        expanded_motif_names.append(f"{motif_name}_alt")
        expanded_motif_names.append(f"{motif_name}_diff")
    columns.extend(expanded_motif_names)

    df = pd.DataFrame(data, columns=columns)
    return df

def read_vcf_to_variant(vcf: pysam.VariantFile) -> Iterator[tuple[Variant, VcfInfo]]:
    for record in vcf:
        chromosome = record.chrom
        position = record.pos
        id_ = "." if not record.id else record.id
        ref = record.ref
        # alt handled below
        quality = "." if not record.qual else record.qual
        filter_ = "." if not (f:="".join(str(f) for f in record.filter)) else f
        info = "." if not (i:=";".join(f"{k}={v}" for k, v in record.info.items())) else i

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
        
        yield (Variant(chromosome, position, "." if ref is None else ref, alt, strand_dir),
            VcfInfo(chromosome, position, id_, ref, alt, quality, filter_, info))

def read_pandas_to_variant(tsv: TextIO) -> Iterator[tuple[Variant, VcfInfo]]:
    headers = tsv.readline().strip().split('\t')  # Adjust index based on required columns
    for line in tsv:
        data = line.strip().split('\t')
        # Basic fields
        chromosome, position, id_, ref, alt, quality, filter_ = data[:7]
        # Combine extra fields into info
        info = ";".join(f"{headers[i]}={data[i]}" for i in range(7, len(headers)))

        # Assuming strand information is available and Variant class is defined properly
        strand = data[7]  # Adjust index for the 'strand' field
        if strand == '+':
            strand_dir = StrandDirection.FORWARD
        elif strand == '-':
            strand_dir = StrandDirection.REVERSE
        else:
            strand_dir = StrandDirection.UNKNOWN

        yield (Variant(chromosome, int(position), ref, alt, strand_dir),
                VcfInfo(chromosome, int(position), id_, ref, alt, quality, filter_, info))

def calculuate_motifs(RBPmotifs: list[RBPsplice], variant_context: VariantContext) -> Optional[list]:
    if (doMSE(variant_context, variant_context.strand_direction == StrandDirection.FORWARD, '5')
        or doMSE(variant_context, variant_context.strand_direction == StrandDirection.FORWARD, '3')):
            return None

    motif_scores = []
    for motif in RBPmotifs:
        # print(motif.name)
        # print(motif.calculate(variant_context.ref_sequence(motif.length)))
        # print(motif.calculate(variant_context.alt_sequence(motif.length)))

        # diff only version
            # a = motif.calculate_variant(variant_context)
            # a = "0" if a == 0 else a
            # motif_scores.append(a)

        # diff and ref version
        ref, alt = motif.calculate_variant_ref_alt(variant_context)
        diff = round(alt - ref, 3)
        ref = round(ref, 3)
        alt = round(alt, 3)

        alt = "0" if alt == 0 else alt
        ref = "0" if ref == 0 else ref
        diff = "0" if diff == 0 else diff

        motif_scores.extend([ref, alt, diff])
    
    return motif_scores

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
    # vcf_variant_iterator = read_vcf_to_variant(pysam.VariantFile(sys.argv[1]))
    tsv_variant_iterator = read_pandas_to_variant(open(sys.argv[1]))
    variant_iterator = tsv_variant_iterator
    output_path = sys.argv[2]
    reference_genome = pysam.FastaFile(sys.argv[3])

    motifs = []
    for _, motif in ESEfinder_motifs.motifs.items():
        motifs.append(RBPsplice.from_2D_list(motif['matrix'], motif['name'], threshold=motif['threshold']))

    for name, motif in RCRUNCH_motifs.motifs.items():
        motifs.append(RBPsplice.from_RCRUNCH(motif, name, threshold=0, arr_by_base=False))

    if not is_path_writable(output_path):
        raise ValueError(f"File path '{output_path}' is not writable!")

    df = calculate_variants(variant_iterator, reference_genome, motifs)

    df.to_csv(output_path, encoding='utf-8', index=False, sep='\t')

if __name__ == "__main__":
    main()
