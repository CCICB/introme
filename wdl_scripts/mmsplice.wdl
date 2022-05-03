version 1.0

task MMSplice {
	input {
		File vcf
		File ref_genome
    File gtf
	}

	command {
		cd ../MMSplice_MTSplice
		python3 run_mmsplice.py --vcf ${vcf} --fasta ${ref_genome} --gtf ${gtf} --output results.mmsplice.vcf
		cp results.mmsplice.vcf /cromwell_root/
	}

	output {
		Array[File] mmsplice_output = glob("*.mmsplice.vcf")
	}

	runtime {
		docker: 'patsul/introme-mmsplice'
		memory: "10G"
		cpu: "4"
	}
}

workflow test {
	input {
		File vcf
		File ref_genome
        File gtf
	}

	call MMSplice {
		input:
			vcf=vcf,
			ref_genome=ref_genome,
			gtf=gtf
	}
}
