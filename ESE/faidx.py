import pysam

def faidx_context(ref_genome: pysam.FastaFile, chrom: str, pos: int,
                  ref_allele: str, context_len: int, prepend_chr=False):
    if prepend_chr and not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    elif not prepend_chr and chrom.startswith("chr"):
        chrom.lstrip("chr")
    
    ref_len = len(ref_allele)

    start = pos - context_len
    end = pos + ref_len + context_len - 1
    sequence = ref_genome.fetch(region=f"{chrom}:{start}-{end}")

    before = sequence[:context_len]
    try:
        assert(ref_allele == sequence[context_len:context_len + ref_len].upper())
    except AssertionError:
        print(f"reference allele ({ref_allele}) did not match the provided reference genome "
              f"({sequence[context_len:context_len + ref_len]}) at {chrom}:{pos}")
        # exit(1)
    after = sequence[context_len + ref_len:]

    return (before, ref_allele, after)

def ref_and_alt(ref_genome: pysam.FastaFile, chrom: str, pos: int, ref_allele: str,
                alt_allele: str, context_len: int, prepend_chr=False):
    before, ref_allele, after = faidx_context(ref_genome, chrom, pos, ref_allele, context_len, prepend_chr)
    # return (f"{before}{ref_allele}{after}", f"{before}{alt_allele}{after}")
    return before, ref_allele, after, alt_allele

def make_window(before: str, ref: str, after: str, alt: str, window_len: int):
    window_len -= 1
    before = before[-window_len:]
    after = after[:window_len]
    # return f"{before}{ref}{after}", f"{before}{alt}{after}"
    return before, ref, after, alt

if __name__ == "__main__":
    import sys
    ref_genome = pysam.FastaFile(sys.argv[1])

    def test_faidx(chrom: str, pos: int, ref_allele: str, alt_allele: str, expected: tuple):
        a = faidx_context(ref_genome, chrom, pos, ref_allele, CONTEXT_LENGTH, prepend_chr=True)
        try:
            assert((a, [len(x) for x in a]) == expected)
        except AssertionError:
            print(f"Assertion that {a}, {[len(x) for x in a]} == {expected} failed")

    def test_ref_and_alt(chrom: str, pos: int, ref_allele: str, alt_allele: str, expected: tuple):
        before, ref_allele, after, alt_allele = ref_and_alt(ref_genome, chrom, pos, ref_allele, alt_allele, CONTEXT_LENGTH, prepend_chr=True)
        a = f"{before}{ref_allele}{after}", f"{before}{alt_allele}{after}"
        try:
            assert(a == expected)
        except AssertionError:
            print(f"Assertion that {a, [len(x) for x in a]} == {expected, [len(x) for x in expected]} "
                    f"failed in case {chrom}:{pos} {ref_allele}>{alt_allele}")        

    # Testing correct splitting of ref allele
    CONTEXT_LENGTH = 7
    test_faidx("1", 77297580, "G", "T", (('GGAAAGG', 'G', 'TACTCAG'), [7, 1, 7]))
    test_faidx("1", 77297580, "G", "GA", (('GGAAAGG', 'G', 'TACTCAG'), [7, 1, 7]))
    test_faidx("1", 77297579, "GG", "TA", (('TGGAAAG', 'GG', 'TACTCAG'), [7, 2, 7]))
    test_faidx("1", 77297580, "GT", "A", (('GGAAAGG', 'GT', 'ACTCAGA'), [7, 2, 7]))

    test_faidx("1", 1393393, "CACCATGAG", "C", (('CTCAACT', 'CACCATGAG', 'GTCTGGA'), [7, 9, 7]))

    test_faidx("1", 42189300, "C", "CA", (('CAACACT', 'C', 'ACCTTGA'), [7, 1, 7]))

    # Testing relevant subsequence for a specific window size can be retrieved
    test_ref_and_alt("1", 77297580, "G", "T",    ('GGAAAGGGTACTCAG',  'GGAAAGGTTACTCAG'))
    test_ref_and_alt("1", 77297580, "G", "GA",   ('GGAAAGGGTACTCAG',  'GGAAAGGGATACTCAG'))
    test_ref_and_alt("1", 77297579, "GG", "TA", ('TGGAAAGGGTACTCAG', 'TGGAAAGTATACTCAG'))
    test_ref_and_alt("1", 77297580, "GT", "A",   ('GGAAAGGGTACTCAGA', 'GGAAAGGAACTCAGA'))
    CONTEXT_LENGTH = 9
    test_ref_and_alt("1", 77297580, "G", "T",    ('GTGGAAAGGGTACTCAGAG',  'GTGGAAAGGTTACTCAGAG'))
    test_ref_and_alt("1", 77297580, "G", "GA",   ('GTGGAAAGGGTACTCAGAG',  'GTGGAAAGGGATACTCAGAG'))
    test_ref_and_alt("1", 77297579, "GG", "TA", ('AGTGGAAAGGGTACTCAGAG', 'AGTGGAAAGTATACTCAGAG'))
    test_ref_and_alt("1", 77297580, "GT", "A",   ('GTGGAAAGGGTACTCAGAGT', 'GTGGAAAGGAACTCAGAGT'))

    res0 = ref_and_alt(ref_genome, "1", 77297580, "G", "T", CONTEXT_LENGTH, prepend_chr=True)
    res1 = ref_and_alt(ref_genome, "1", 77297579, "GG", "TA", CONTEXT_LENGTH, prepend_chr=True)
    print(res0)

    seq1, seq2 = "".join(res0[0:3]), "".join(res0[0] + res0[3] + res0[2])
    print(seq1, seq2)
    # Test for win8
    window = 8
    for i in range(window):
        offset = i + (CONTEXT_LENGTH - window) + 1
        print(offset, seq1[offset:offset + window])
        print(offset, seq2[offset:offset + window])
    x = (CONTEXT_LENGTH - window) + 1
    print(seq1[x:x + (2 * window - 1)])
    print(seq2[x:x + (2 * window - 1)])    

    print("=====")

    print(make_window(*res0, 5))
    print(make_window(*res1, 5))
    print(make_window(*res0, 8))
    print(make_window(*res1, 8))

    # window = 7
    # for i in range(window):
    #     offset = i + (CONTEXT_LENGTH - window) + 1
    #     print(offset, seq1[offset:offset + window])    
    # x = (CONTEXT_LENGTH - window) + 1
    # print(seq1[x:x + (2 * window - 1)])    

    # print("=====")

    # window = 6
    # for i in range(window):
    #     offset = i + (CONTEXT_LENGTH - window) + 1
    #     print(offset, seq1[offset:offset + window])    
    # x = (CONTEXT_LENGTH - window) + 1
    # print(seq1[x:x + (2 * window - 1)])    
