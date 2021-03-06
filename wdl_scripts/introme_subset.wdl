version 1.0

task subset {
	input {
		File vcf
		File bam
		String prefix
	}

	command {
		./subsetting.sh -v ${vcf} -b ${bam} -p ${prefix}
	}

	output {
		Array[File] subsetting_output = glob("*.vcf")
	}

	runtime {
		docker: 'patsul/introme-subset'
	}
}

workflow introme_subset {
	input {
		File vcf
		File bam
	}

	call vcfanno {
		input:
			vcf=vcf,
			bam=bam,
			prefix=prefix
	}
}
