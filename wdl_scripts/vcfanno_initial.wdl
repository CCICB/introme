version 1.0

task vcfanno {
	input {
		File vcf
		File gnomad
		File gnomad_tbi
		File cadd
		File cadd_tbi
		File MGRB
		File MGRB_tbi
		File dbscSNV
		File dbscSNV_tbi
		File SPIDEX
		File SPIDEX_tbi
		File branchpointer
		File branchpointer_tbi
		File gencode
		File gencode_tbi
		File regions
		File regions_tbi
		File locations
		File locations_tbi
		File U12
		File U12_tbi
	}

	command {
		cd /app		
		vcfanno -base-path $(echo ${regions} | cut -f 1-$(echo ${regions} | grep -o '/' | grep -c .) -d '/') conf_test.toml ${vcf} > /cromwell_root/vcfanno_test_output.vcf
	}

	output {
		Array[File] vcfanno_output = glob("*.vcf")
	}

	runtime {
		docker: 'patsul/introme-vcfanno'
	}
}

workflow test {
	input {
		File vcf
		File gnomad
		File gnomad_tbi
		File cadd
		File cadd_tbi
		File MGRB
		File MGRB_tbi
		File dbscSNV
		File dbscSNV_tbi
		File SPIDEX
		File SPIDEX_tbi
		File branchpointer
		File branchpointer_tbi
		File gencode
		File gencode_tbi
		File regions
		File regions_tbi
		File locations
		File locations_tbi
		File U12
		File U12_tbi
	}

	call vcfanno {
		input:
			vcf=vcf,
			gnomad=gnomad,
			gnomad_tbi=gnomad_tbi,
			cadd=cadd,
			cadd_tbi=cadd_tbi,
			MGRB=MGRB,
			MGRB_tbi=MGRB_tbi,
			dbscSNV=dbscSNV,
			dbscSNV_tbi=dbscSNV_tbi,
			SPIDEX=SPIDEX,
			SPIDEX_tbi=SPIDEX_tbi,
			branchpointer=branchpointer,
			branchpointer_tbi=branchpointer_tbi,
			gencode=gencode,
			gencode_tbi=gencode_tbi,
			regions=regions,
			regions_tbi=regions_tbi,
			locations=locations,
			locations_tbi=locations_tbi,
			U12=U12,
			U12_tbi=U12_tbi
	}
}
