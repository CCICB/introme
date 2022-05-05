version 1.0

task Spliceogen {
	input {
		File vcf
		File ref_genome
		File gtf
	}

	command {
		cd ../Spliceogen
		./RUN.sh -input ${vcf} -fasta ${ref_genome} -gtf ${gtf}
		cp output/*.vcf_out.txt /cromwell_root/
	}
	
    output {
		Array[File] spliceogen_output = glob("*.vcf_out.txt")
    }
    
    runtime {
    	docker: 'patsul/introme-spliceogen'
    }
}

workflow run_spliceogen {
	input {
    		File vcf
        	File ref_genome
        	File gtf
    }
    
    call Spliceogen {
		input:
			vcf=vcf,
			ref_genome=ref_genome,
			gtf=gtf
    }
}
