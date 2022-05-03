version 1.0

task filter {
	input {
		File annotated_vcf
	}

	command {
		./filter.sh -v ${vcf}
	}

	output {
		Array[File] filtered_output = glob("*.vcf")
	}

	runtime {
		docker: 'patsul/introme-filter'
	}
}

workflow introme_filter {
	input {
		File annotated_vcf
	}

	call filter {
		input:
			annotated_vcf=annotated_vcf,
	}
}
