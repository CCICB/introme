version 1.0

task SpliceAI {
	input {
		File vcf
		File ref_genome
		Int grchX
		Int distance
		Int mask
	}

	command {
		spliceai -I ${vcf} -O spliceai.vcf -R ${ref_genome} -A grch${grchX} -D ${distance} -M ${mask}
	}

	output {
		Array[File] spliceai_output = glob("*.vcf")
	}

	runtime {
		docker: 'patsul/introme-spliceai'
	}
}

workflow run_spliceai {
	input {
		File vcf
		File ref_genome
		Int grchX
		Int distance
		Int mask
	}

	call SpliceAI {
		input:
			vcf=vcf,
			ref_genome=ref_genome,
			grchX=grchX,
			distance=distance,
			mask=mask
	}
}
