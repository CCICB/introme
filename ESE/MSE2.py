from varconv.variants import Variant, VariantContext
from typing import Literal
import timeit
import pysam
import sys

from maxentpy import maxent
matrix3 = maxent.load_matrix3()
matrix5 = maxent.load_matrix5()

variants = [
    Variant("chr1", 77297580, "G", "T"),
]

contexts = [VariantContext(pysam.FastaFile(sys.argv[1]), var, 22) for var in variants]

def sliding_window(string, n):
    # Ensure the window size is valid
    if n > len(string) or n <= 0:
        raise ValueError("Window size must be positive and less than or equal to the string length")

    # Get the number of windows
    num_windows = len(string) - n + 1

    # Generate and return the windows
    return [string[i:i+n] for i in range(num_windows)]


THRESHOLD = 0

def doMSE(variant_context: VariantContext, forward_strand: bool, splicing_end: Literal['3', '5'], context_length: int):
    ref_max = -999
    alt_at_ref_max = -999
    ref_max_pos = -999
    alt_max = -999
    ref_at_alt_max = -999
    alt_max_pos = -999

    if (splicing_end == '3'):
        context_length = 23
        fn = maxent.score3
        matrix = matrix3
    elif (splicing_end == '5'):
        context_length = 9
        fn = maxent.score5
        matrix = matrix5
    else:
        raise ValueError()

    if (forward_strand):
        print(f"SLIDING WIN {context_length}")
        gen = enumerate(
            zip(sliding_window(variant_context.ref_sequence_fwd(context_length), context_length),
                sliding_window(variant_context.alt_sequence_fwd(context_length), context_length))
        )
    else:
        gen = enumerate(
                zip(sliding_window(variant_context.ref_sequence_rev(context_length), context_length),
                    sliding_window(variant_context.alt_sequence_rev(context_length), context_length))
            )

    for i, (ref, alt) in gen:
        ref_score = fn(ref, matrix)
        alt_score = fn(alt, matrix)

        if ref_score > ref_max:
            ref_max = ref_score
            alt_at_ref_max = alt_score
            ref_max_pos = i

        if alt_score > alt_max:
            alt_max = alt_score
            ref_at_alt_max = ref_score
            alt_max_pos = i

        print(f"{ref}\n{alt}, {alt_score:.2f} - {ref_score:.2f} = {alt_score - ref_score:.2f}")
        if ref_score > THRESHOLD or alt_score > THRESHOLD:
            print("", ref, "\n", alt, ref_score, alt_score)

    print(f'{ref_max=:.2f} {alt_at_ref_max=:.2f} at i={ref_max_pos}  {alt_max=:.2f} {ref_at_alt_max=:.2f} at i={alt_max_pos}')

for c in contexts:
    doMSE(c, True, '3', 23)
    doMSE(c, True, '5', 9)