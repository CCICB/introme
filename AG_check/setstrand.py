#!/usr/bin/env python3

import pysam
import random

input_vcf = "AG_check/files/pedcbioportal_short_500.vcf"
output_vcf = "AG_check/files/pedcbioportal_strand_500.vcf"

vcf_in = pysam.VariantFile(input_vcf)
vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

for record in vcf_in.fetch():
    random_number1 = random.random()
    random_number2 = random.random()

    strands = []
    if random_number1 <= 0.25:
        strands.append("+1")
    if random_number2 <= 0.25:
        strands.append("-1")

    if strands:
        record.info["strand"] = ",".join(strands)

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
