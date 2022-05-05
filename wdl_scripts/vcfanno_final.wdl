version 1.0

task vcfanno {
	input {
		File vcf
		File spliceai
		File spliceai_tbi
		File mmsplice
		File mmsplice_tbi
		File spliceogen
		File spliceogen_tbi
	}

	command {
		cd /app		
		vcfanno -base-path $(echo ${spliceai} | cut -f 1-$(echo ${spliceai} | grep -o '/' | grep -c .) -d '/') conf_test.toml ${vcf} > /cromwell_root/vcfanno_test_output.vcf
	}

	output {
		Array[File] vcfanno_output = glob("*.vcf")
	}

	runtime {
		docker: 'patsul/introme-vcfanno'
	}
}

workflow run_vcfanno {
	input {
		File vcf
		File spliceai
		File spliceai_tbi
		File mmsplice
		File mmsplice_tbi
		File spliceogen
		File spliceogen_tbi
	}

	call vcfanno {
		input:
			vcf=vcf,
			spliceai=spliceai,
			spliceai_tbi=spliceai_tbi,
			mmsplice=mmsplice,
			mmsplice_tbi=mmsplice_tbi,
			spliceogen=spliceogen,
			spliceogen_tbi=spliceogen_tbi
	}
}
