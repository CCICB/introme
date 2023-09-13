import pysam
from dataclasses import dataclass
from Bio.Seq import reverse_complement

# Variants with ref/alt allele and context window before and after
@dataclass(frozen=True)
class VariantContext():
    before: str
    ref_allele: str
    alt_allele: str
    after: str

    context_length: int

    # TODO: find a way to not have to remake strings for context window
    # TODO: assert context_length is long enough for window_length
    def ref_sequence(self, window_length: int=None) -> str:
        if window_length is not None:
            window_length -= 1
            return f"{self.before[-window_length:]}{self.ref_allele}{self.after[:window_length]}"
        return f"{self.before}{self.ref_allele}{self.after}"
    
    def alt_sequence(self, window_length: int=None) -> str:
        if window_length is not None:
            window_length -= 1
            return f"{self.before[-window_length:]}{self.alt_allele}{self.after[:window_length]}"
        return f"{self.before}{self.alt_allele}{self.after}"
    
    def ref_sequence_rev(self, window_length: int=None) -> str:
        return reverse_complement(self.ref_sequence(window_length))

    def alt_sequence_rev(self, window_length: int=None) -> str:
        return reverse_complement(self.alt_sequence(window_length))

    def debug(self) -> tuple:
        tup = (self.before, self.ref_allele, self.after, self.alt_allele)
        return tup, [len(elem) for elem in tup]
    
    def __repr__(self) -> str:
        return str(self.debug())

# Variants extracted from vcf file
@dataclass(frozen=True)
class Variant():
    chrom: str
    position: int
    ref_allele: str
    alt_allele: str
    
    def faidx_context(self, ref_genome: pysam.FastaFile, context_len: int) -> VariantContext:
        ref_len = len(self.ref_allele)
        start = self.position - context_len
        end = self.position + ref_len + context_len - 1
        sequence = ref_genome.fetch(region=f"{self.chrom}:{start}-{end}").upper()

        before = sequence[:context_len]
        try:
            assert(self.ref_allele == sequence[context_len:context_len + ref_len])
        except AssertionError:
            print(f"Variant reference allele ({self.ref_allele}) did not match the provided reference genome "
                f"({sequence[context_len:context_len + ref_len]}) at {self.chrom}:{self.position}")
            # exit(1)
            # pass
        after = sequence[context_len + ref_len:]

        return VariantContext(before, self.ref_allele, self.alt_allele, after, context_len)

def main():
    pass

if __name__ == '__main__':
    main()
