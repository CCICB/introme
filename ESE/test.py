from introme.ESE.variants import Variant, VariantContext
import pysam

def test_faidx(ref_genome: pysam.FastaFile):
    def test(variant: Variant, context_length: int, expected: tuple):
        variant_context = variant.faidx_context(ref_genome, context_length)
        try:
            assert(variant_context.debug() == expected)
        except AssertionError:
            print(f"Assertion that {variant_context} == {expected} failed")

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

def main():
    import sys
    ref_genome = pysam.FastaFile(sys.argv[1])
    test_faidx(ref_genome)

if __name__ == '__main__':
    main()