from variants import Variant, VariantContext
import pysam
import Bio.Seq

def test_faidx(ref_genome: pysam.FastaFile):
    def test(variant: Variant, context_length: int, expected: tuple):
        variant_context = variant.faidx_context(ref_genome, context_length)
        try:
            assert(variant_context.debug() == expected)
        except AssertionError:
            print(f"Assertion that {variant_context} == {expected} failed  for {variant}")

    CONTEXT_LENGTH = 7
    test(Variant("chr1", 77297580, "G", "T"), CONTEXT_LENGTH, (('GGAAAGG', 'G', 'TACTCAG', 'T'), [7, 1, 7, 1]))
    test(Variant("chr1", 77297580, "G", "GA"), CONTEXT_LENGTH, (('GGAAAGG', 'G', 'TACTCAG', 'GA'), [7, 1, 7, 2]))
    test(Variant("chr1", 77297579, "GG", "TA"), CONTEXT_LENGTH, (('TGGAAAG', 'GG', 'TACTCAG', 'TA'), [7, 2, 7, 2]))
    test(Variant("chr1", 77297580, "GT", "A"), CONTEXT_LENGTH, (('GGAAAGG', 'GT', 'ACTCAGA', 'A'), [7, 2, 7, 1]))

    test(Variant("chr1", 1393393, "CACCATGAG", "C"), CONTEXT_LENGTH, (('CTCAACT', 'CACCATGAG', 'GTCTGGA', 'C'), [7, 9, 7, 1]))
    test(Variant("chr1", 42189300, "C", "CA"), CONTEXT_LENGTH, (('CAACACT', 'C', 'ACCTTGA', 'CA'), [7, 1, 7, 2]))

    CONTEXT_LENGTH = 9
    test(Variant("chr1", 77297580, "G", "T"), CONTEXT_LENGTH, (('GTGGAAAGG', 'G', 'TACTCAGAG', 'T'),  [9, 1, 9, 1]))
    test(Variant("chr1", 77297580, "G", "GA"), CONTEXT_LENGTH, (('GTGGAAAGG', 'G', 'TACTCAGAG', 'GA'),  [9, 1, 9, 2]))
    test(Variant("chr1", 77297579, "GG", "TA"), CONTEXT_LENGTH, (('AGTGGAAAG', 'GG', 'TACTCAGAG', 'TA'), [9, 2, 9, 2]))
    test(Variant("chr1", 77297580, "GT", "A"), CONTEXT_LENGTH, (('GTGGAAAGG', 'GT', 'ACTCAGAGT', 'A'), [9, 2, 9, 1]))

def test_ref_alt_sequence(ref_genome: pysam.FastaFile):
    from dataclasses import dataclass
    @dataclass
    class Expected(dict):
        ref_seq_fwd: str
        alt_seq_fwd: str

    def test_forward(variant: Variant, context_length: int, expected: Expected):
        variant_context = variant.faidx_context(ref_genome, context_length)
        try:
            assert(variant_context.ref_sequence_fwd() == expected.ref_seq_fwd)
        except AssertionError:
            print(f"Assertion that {variant_context.ref_sequence_fwd()} == {expected.ref_seq_fwd} failed for {variant}")

        try:
            assert(variant_context.alt_sequence_fwd() == expected.alt_seq_fwd)
        except AssertionError:
            print(f"Assertion that {variant_context.alt_sequence_fwd()} == {expected.alt_seq_fwd} failed for {variant}")

    def test_reverse(variant: Variant, context_length: int, expected: Expected):
        variant_context = variant.faidx_context(ref_genome, context_length)
        try:
            assert(variant_context.ref_sequence_rev() == Bio.Seq.reverse_complement(expected.ref_seq_fwd))
        except AssertionError:
            print(f"Assertion that {variant_context.ref_sequence_rev()} == "
                  f"{Bio.Seq.reverse_complement(expected.ref_seq_fwd)} failed for {variant}")

        try:
            assert(variant_context.alt_sequence_rev() == Bio.Seq.reverse_complement(expected.alt_seq_fwd))
        except AssertionError:
            print(f"Assertion that {variant_context.alt_sequence_rev()} == "
                  f"{Bio.Seq.reverse_complement(expected.alt_seq_fwd)} failed for {variant}")

    # Forward strand
    CONTEXT_LENGTH = 7
    test_forward(Variant("chr1", 77297580, "G", "T"), CONTEXT_LENGTH,
             Expected('GGAAAGG' + 'G' + 'TACTCAG', 'GGAAAGG' + 'T' + 'TACTCAG'))
    test_forward(Variant("chr1", 77297580, "G", "GA"), CONTEXT_LENGTH,
             Expected('GGAAAGG' + 'G' + 'TACTCAG','GGAAAGG' + 'GA' + 'TACTCAG'))
    test_forward(Variant("chr1", 77297579, "GG", "TA"), CONTEXT_LENGTH,
             Expected('TGGAAAG' + 'GG' + 'TACTCAG','TGGAAAG' + 'TA' + 'TACTCAG'))
    test_forward(Variant("chr1", 77297580, "GT", "A"), CONTEXT_LENGTH,
             Expected('GGAAAGG' + 'GT' + 'ACTCAGA','GGAAAGG' + 'A' + 'ACTCAGA'))

    # Reverse strand, supply Expected() as forward strand sequence, sourced from ucsc DNA browser
    test_reverse(Variant("chr1", 5873238, "C", "T"), CONTEXT_LENGTH,
             Expected('TGTCCTC' + 'C' + 'GTTGCCC', 'TGTCCTC' + 'T' + 'GTTGCCC'))
    test_reverse(Variant("chr1", 1393393, "CACCATGAG", "C"), CONTEXT_LENGTH,
             Expected('CTCAACT' + 'CACCATGAG' + 'GTCTGGA', 'CTCAACT' + 'C' + 'GTCTGGA'))
    test_reverse(Variant("chr1", 42189300, "C", "CA"), CONTEXT_LENGTH,
             Expected('CAACACT' + 'C' + 'ACCTTGA', 'CAACACT' + 'CA' + 'ACCTTGA'))



def main():
    import sys
    ref_genome = pysam.FastaFile(sys.argv[1])
    test_faidx(ref_genome)
    test_ref_alt_sequence(ref_genome)

if __name__ == '__main__':
    main()
