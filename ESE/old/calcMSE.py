import pysam
import sys
from varconv.variants import Variant, VariantContext, StrandDirection
# from variants import Variant, VariantContext, StrandDirection
from scoring import read_pandas_to_variant, VcfInfo
from typing import Callable, Iterator
from MSE import doMSE

CONTEXT_LENGTH = 22

def add_hello_world_decorator(original_function: Callable, fasta: pysam.FastaFile) -> Callable:
    def report_with_prefix(pre: str, vcf_info: VcfInfo, res: list):
        vcf_info.addInfo(f"{pre}effect", res[0])

        vcf_info.addInfo(f"{pre}ref_max", round(res[1], 2))
        vcf_info.addInfo(f"{pre}alt_at_ref_max", round(res[2], 2))
        vcf_info.addInfo(f"{pre}ref_max_pos", round(res[3], 2))

        vcf_info.addInfo(f"{pre}alt_max", round(res[4], 2))
        vcf_info.addInfo(f"{pre}ref_at_alt_max", round(res[5], 2))
        vcf_info.addInfo(f"{pre}alt_max_pos", round(res[6], 2))
        
    def decorated_function(*args, **kwargs) -> Iterator[tuple[Variant, VcfInfo]]:
        for variant, vcf_info in original_function(*args, **kwargs):
            context = VariantContext(fasta, variant, CONTEXT_LENGTH)
            # context = variant.faidx_context(fasta, CONTEXT_LENGTH)
            # print(context.variant.strand_direction, context.variant.strand_direction == StrandDirection.FORWARD)
            res3 = doMSE(context, context.variant.strand_direction == StrandDirection.FORWARD, '3')
            res5 = doMSE(context, context.variant.strand_direction == StrandDirection.FORWARD, '5')
            # res3 = doMSE(context, True, '3')
            # res5 = doMSE(context, True, '5')

            report_with_prefix('3_', vcf_info, res3)
            report_with_prefix('5_', vcf_info, res5)
            yield variant, vcf_info
    return decorated_function

if __name__ == "__main__":
    fasta = pysam.FastaFile(sys.argv[1])

    read_pandas_to_variant_decorated = add_hello_world_decorator(read_pandas_to_variant, fasta)
    variant_iterator = read_pandas_to_variant_decorated(open(sys.argv[2]))

    for variant, vcf_info in variant_iterator:
        # variant and vcf_info will have the additional "hello=world" in the info
        print(vcf_info.toString())