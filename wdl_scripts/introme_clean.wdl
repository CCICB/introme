version 1.0

task clean {
	input {
		File annotated_vcf
	}

	command {
		./clean.sh -v ${vcf}
	}

	output {
		Array[File] cleaned_output = glob("*.vcf")
	}

	runtime {
		docker: 'patsul/introme-clean'
	}
}

workflow introme_clean {
	input {
		File annotated_vcf
	}

	call clean {
		input:
			annotated_vcf=annotated_vcf,
	}
}
