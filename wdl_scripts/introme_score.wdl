version 1.0

task score {
	input {
		File final_annotated_vcf
	}

	command {
		./score.sh -v ${final_annotations_vcf}
	}

	output {
		Array[File] introme_output = glob("*.tsv")
	}

	runtime {
		docker: 'patsul/introme-score'
	}
}

workflow introme_score {
	input {
		File final_annotated_vcf
	}

	call score {
		input:
			final_annotated_vcf=final_annotated_vcf,
	}
}
