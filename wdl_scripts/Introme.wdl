version 1.0

## WORKFLOW DEFINITIONS 
workflow introme {
	input {
		File vcf
		File ref_genome
		File gtf
		Int grchX
		File bed
		String prefix
		Int quality_filter
		String prefix
		File gnomad
		File gnomad_tbi
		File MGRB
		File MGRB_tbi
		File regions
		File regions_tbi
		Float allele_frequency
		File cadd
		File cadd_tbi
		File dbscSNV
		File dbscSNV_tbi
		File SPIDEX
		File SPIDEX_tbi
		File branchpointer
		File branchpointer_tbi
		File gencode
		File gencode_tbi
		File locations
		File locations_tbi
		File U12
		File U12_tbi
	}

	call subset {
		input:
			vcf=vcf,
			bed=bed,
			prefix=prefix,
			quality_filter=quality_filter
	}

	call vcfanno_freq {
		input:
			prefix=prefix,
			vcf=subset.out_filter,
			gnomad=gnomad,
			gnomad_tbi=gnomad_tbi,
			MGRB=MGRB,
			MGRB_tbi=MGRB_tbi,
			regions=regions,
			regions_tbi=regions_tbi,
			gencode=gencode,
			gencode_tbi=gencode_tbi
	}
    
	call filter {
		input:
			annotated_vcf=vcfanno_freq.vcfanno_output,
			prefix=prefix,
			allele_frequency=allele_frequency
	}

	call SpliceAI {
		input:
			vcf=filter.out_rmanno,
			ref_genome=ref_genome,
			grchX=grchX,
			distance=1000,
			mask=0,
			prefix=prefix
	}
    
	call MMSplice {
		input:
			vcf=filter.out_rmanno,
			ref_genome=ref_genome,
			gtf=gtf,
			prefix=prefix
	}
    
	call Spliceogen {
		input:
			vcf=filter.out_rmanno,
			ref_genome=ref_genome,
			gtf=gtf,
			prefix=prefix
    }

	call vcfanno_splicing_anno {
		input:
			vcf=filter.out_filter,
			prefix=prefix,
			cadd=cadd,
			cadd_tbi=cadd_tbi,
			dbscSNV=dbscSNV,
			dbscSNV_tbi=dbscSNV_tbi,
			SPIDEX=SPIDEX,
			SPIDEX_tbi=SPIDEX_tbi,
			branchpointer=branchpointer,
			branchpointer_tbi=branchpointer_tbi,
			locations=locations,
			locations_tbi=locations_tbi,
			U12=U12,
			U12_tbi=U12_tbi
	}
    
	call extract {
		input:
			annotated_vcf=filter.out_filter,
			prefix=prefix,
			spliceai_output=SpliceAI.spliceai_output,
			mmsplice_output=MMSplice.mmsplice_output,
			spliceogen_output=Spliceogen.spliceogen_output,
			ref_genome=ref_genome
	}
    
	call vcfanno_splicing {
		input:
			vcf=vcfanno_splicing_anno.vcfanno_spliceanno,
			prefix=prefix,
			spliceai=extract.out_spliceai,
			spliceai_tbi=extract.out_spliceai_tbi,
			spliceogen=extract.out_spliceogen,
			spliceogen_tbi=extract.out_spliceogen_tbi,
			mmsplice=extract.out_mmsplice,
			mmsplice_tbi=extract.out_mmsplice_tbi,
			functions=extract.out_functions,
			functions_tbi=extract.out_functions_tbi
	}
    
	call clean {
		input:
			vcf=vcfanno_splicing.vcfanno_splicing,
			prefix=prefix
	}

	call score {
		input:
			tsv=clean.clean_output,
			prefix=prefix
	}
}


## TASK DEFINITIONS 
task subset {
	input {
		File vcf
		File bed
		String prefix
		Int? quality_filter = 1
	}

	command {
		cd ..
		./subsetting.sh -v ${vcf} -b ${bed} -p ${prefix} -q ${quality_filter}
		cp *.vcf.gz /cromwell_root/
	}

	output {
		File out_subset = "${prefix}.subset.vcf.gz"
		File out_filter = "${prefix}.subset.highquality.vcf.gz"
	}

	runtime {
		docker: 'patsul/introme-base'
	}
}


task vcfanno_freq {
	input {
		String prefix
		File vcf
		File gnomad
		File gnomad_tbi
		File MGRB
		File MGRB_tbi
		File regions
		File regions_tbi
		File gencode
		File gencode_tbi
	}

	command {
		cd /app		
		vcfanno -base-path $(echo ${gnomad} | cut -f 1-$(echo ${gnomad} | grep -o '/' | grep -c .) -d '/') conf_frequency.toml ${vcf} | bgzip > /cromwell_root/${prefix}.subset.highquality.annotated.vcf.gz
	}

	output {
		File vcfanno_output = "${prefix}.subset.highquality.annotated.vcf.gz"
	}

	runtime {
		docker: 'patsul/introme-vcfanno'
	}
}


task filter {
	input {
		File annotated_vcf
		String prefix
		Float? allele_frequency = 0.01
	}

	command {
		cd ..
		./filter.sh -v ${annotated_vcf} -p ${prefix} -a ${allele_frequency}
		cp *.vcf.gz /cromwell_root/
	}

	output {
		File out_filter = "${prefix}.subset.highquality.annotated.filtered.vcf.gz"
		File out_rmanno = "${prefix}.subset.highquality.annotated.filtered.rmanno.vcf.gz"
	}

	runtime {
		docker: 'patsul/introme-base'
	}
}

task SpliceAI {
	input {
		File vcf
		File ref_genome
		Int grchX
		Int? distance = 1000
		Int? mask = 0
		String prefix
	}

	command {
		gunzip -c ${vcf} > ${prefix}.spliceai_input.vcf
		spliceai -I ${prefix}.spliceai_input.vcf -O /cromwell_root/${prefix}.spliceai.vcf -R ${ref_genome} -A grch${grchX} -D ${distance} -M ${mask}
	}

	output {
		File spliceai_output = "${prefix}.spliceai.vcf"
	}

	runtime {
		docker: 'patsul/introme-spliceai'
		cpu: "32"
	}
}

task MMSplice {
	input {
		File vcf
		File ref_genome
		File gtf
		String prefix
	}

	command {
		cd ../MMSplice_MTSplice
		python3 run_mmsplice.py --vcf ${vcf} --fasta ${ref_genome} --gtf ${gtf} --output ${prefix}.mmsplice.vcf
		cp ${prefix}.mmsplice.vcf /cromwell_root/
	}

	output {
		File mmsplice_output = "${prefix}.mmsplice.vcf"
	}

	runtime {
		docker: 'patsul/introme-mmsplice'
		memory: "10G"
		cpu: "4"
	}
}

task Spliceogen {
	input {
		File vcf
		File ref_genome
		File gtf
		String prefix
	}

	command {
		cd ../Spliceogen
		gunzip -c ${vcf} > ${prefix}.spliceogen_input.vcf
		./RUN.sh -input ${prefix}.spliceogen_input.vcf -fasta ${ref_genome} -gtf ${gtf}
		cp output/*.vcf_out.txt /cromwell_root/${prefix}.spliceogen.txt
	}
	
	output {
		File spliceogen_output = "${prefix}.spliceogen.txt"
	}
    
	runtime {
    		docker: 'patsul/introme-spliceogen'
    }
}

task extract {
	input {
		File spliceai_output
		File mmsplice_output
		File spliceogen_output
		File ref_genome
		File annotated_vcf
		String prefix
	}

	command {
		cd ..
		./extract.sh -a ${spliceai_output} -b ${spliceogen_output} -c ${mmsplice_output} -p ${prefix}
		./AG_GT_check.sh -r ${ref_genome} -v ${annotated_vcf}
		cp introme_annotate.*.vcf.gz /cromwell_root/
		cp introme_annotate.*.vcf.gz.tbi /cromwell_root/
	}

	output {
		File out_spliceogen = "introme_annotate.spliceogen.vcf.gz"
		File out_spliceogen_tbi = "introme_annotate.spliceogen.vcf.gz.tbi"
		File out_spliceai = "introme_annotate.spliceai.vcf.gz"
		File out_spliceai_tbi = "introme_annotate.spliceai.vcf.gz.tbi"
		File out_mmsplice = "introme_annotate.mmsplice.vcf.gz"
		File out_mmsplice_tbi = "introme_annotate.mmsplice.vcf.gz.tbi"
		File out_functions = "introme_annotate.functions.vcf.gz"
		File out_functions_tbi = "introme_annotate.functions.vcf.gz.tbi"
	}

	runtime {
		docker: 'patsul/introme-base'
	}
}

task vcfanno_splicing_anno {
	input {
		File vcf
		String prefix
		File cadd
		File cadd_tbi
		File dbscSNV
		File dbscSNV_tbi
		File SPIDEX
		File SPIDEX_tbi
		File branchpointer
		File branchpointer_tbi
		File locations
		File locations_tbi
		File U12
		File U12_tbi
	}

	command {
		cd /app		
		vcfanno -base-path $(echo ${U12} | cut -f 1-$(echo ${U12} | grep -o '/' | grep -c .) -d '/') conf_splicing_anno.toml ${vcf} | bgzip > /cromwell_root/${prefix}.subset.highquality.annotated.filtered.scored_anno.vcf.gz
	}

	output {
		File vcfanno_spliceanno = "${prefix}.subset.highquality.annotated.filtered.scored_anno.vcf.gz"
	}

	runtime {
		docker: 'patsul/introme-vcfanno'
		disks: "local-disk 75 SSD"
	}
}


task vcfanno_splicing {
	input {
		File vcf
		String prefix
		File spliceai
		File spliceai_tbi
		File spliceogen
		File spliceogen_tbi
		File mmsplice
		File mmsplice_tbi
		File functions
		File functions_tbi
	}

	command {
		cd /app		
		vcfanno -base-path $(echo ${spliceai} | cut -f 1-$(echo ${spliceai} | grep -o '/' | grep -c .) -d '/') conf_splicing_run.toml ${vcf} | bgzip > /cromwell_root/${prefix}.subset.highquality.annotated.filtered.scored.vcf.gz
	}

	output {
		File vcfanno_splicing = "${prefix}.subset.highquality.annotated.filtered.scored.vcf.gz"
	}

	runtime {
		docker: 'patsul/introme-vcfanno'
	}
}

task clean {
	input {
		File vcf
		String prefix
	}

	command {
		cd ..
		./clean.sh -p ${prefix} -v ${vcf}
		cp ${prefix}.annotated.cleaned.tsv /cromwell_root/${prefix}.annotated.cleaned.tsv
	}
	
    output {
		File clean_output = "${prefix}.annotated.cleaned.tsv"
    }
    
    runtime {
		docker: 'patsul/introme-base'
    }
}

task score {
	input {
		File tsv
		String prefix
	}

	command {
		cd ..
		Rscript --vanilla consensus_scoring.R full ${tsv} /cromwell_root/${prefix}.introme.predictions.tsv
	}

	output {
		File out_score = "${prefix}.introme.predictions.tsv"
	}

	runtime {
		docker: 'patsul/introme-scoring'
		memory: "10G"
	}
}
