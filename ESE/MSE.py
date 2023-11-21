from varconv.variants import Variant, VariantContext
import timeit
import pysam
import sys

# from maxentpy import maxent
# from maxentpy.maxent import load_matrix5, load_matrix3, maxent_fast


from maxentpy import maxent  # use normal version of maxent
print(maxent.score5('cagGTAAGT'))
print(maxent.score3('ttccaaacgaacttttgtAGgga'))  # 20 bases in the intron and 3 base in the exon


from maxentpy.maxent import load_matrix5, load_matrix3
def x():
    matrix3 = load_matrix3()
    maxent.score3('ttccaaacgaacttttgtAGgga', matrix=matrix3)

print(timeit.timeit(lambda: x(), number=1))

matrix3 = load_matrix3()
matrix5 = load_matrix5()
print(timeit.timeit(lambda: maxent.score3('ttccaaacgaacttttgtAGgga', matrix=matrix3), number=1000))


variants = [
    Variant("chr1", 77297580, "G", "T"),
    # Variant("chr1", 5873238, "C", "T"),
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


THRESHOLD = 3

for variant_context in contexts:
    ref_max = -999
    alt_at_ref_max = -999
    ref_max_pos = -999
    alt_max = -999
    ref_at_alt_max = -999
    alt_max_pos = -999

    print("3'", variant_context.alt_sequence_fwd(), variant_context.ref_sequence_fwd())
    for i, (ref, alt) in enumerate(
                        zip(sliding_window(variant_context.ref_sequence_fwd(), 23),
                            sliding_window(variant_context.alt_sequence_fwd(), 23))
    ):
        ref_score = maxent.score3(ref, matrix=matrix3)
        alt_score = maxent.score3(alt, matrix=matrix3)

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
            print ("", ref, "\n", alt, ref_score, alt_score)

    print(f'{ref_max=:.2f} {alt_at_ref_max=:.2f} at i={ref_max_pos}  {alt_max=:.2f} {ref_at_alt_max=:.2f} at i={alt_max_pos}')

    ref_max = -999
    alt_at_ref_max = -999
    ref_max_pos = -999
    alt_max = -999
    ref_at_alt_max = -999
    alt_max_pos = -999

    print("5'")
    for ref, alt in zip(sliding_window(variant_context.ref_sequence_fwd(9), 9), sliding_window(variant_context.alt_sequence_fwd(9), 9)):
        ref_score = maxent.score5(ref, matrix=matrix5)
        alt_score = maxent.score5(alt, matrix=matrix5)
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
            print ("", ref, "\n", alt, ref_score, alt_score)

    print(f'{ref_max=:.2f} {alt_at_ref_max=:.2f} at i={ref_max_pos}  {alt_max=:.2f} {ref_at_alt_max=:.2f} at i={alt_max_pos}')



