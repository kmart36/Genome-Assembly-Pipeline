#!/usr/bin/env nextflow

pb_files = Channel.fromPath("/home/kam071/SampleReads/Nymphalis_polychloros/pac_bio/*.fastq", type: 'file')

hc_files = Channel.fromFilePairs("/home/kam071/SampleReads/Nymphalis_polychloros/hi-c/*_{1,2}.fastq", checkIfExists:true)

pb_files.into{files_fastQC; files_jellyfish; files_hifiasm}

hc_files.into{files_hc_hifiasm}

process fastQC {
	conda 'fastqc'

	publishDir '/SampleReads/Nymphalis_polychloros'	

	input:
	file sample from files_fastQC
	
	script:
	"""
		mkdir -p ${sample.baseName}_fastqc
		fastqc --outdir ${sample.baseName}_fastqc \
		${sample}
	"""
}

process jellyfish {
	
	conda 'jellyfish'
	memory '128 GB'
	executor 'slurm'

        input:
	file sample from files_jellyfish
	
	output:
        file "${sample.baseName}.jf" into jellyfiles

	script:
	"""
		jellyfish count -C -m 21 -s 16000000000 -t 10 ${sample} -o ${sample.baseName}.jf

	"""
}

process genomescope {
	
	conda 'jellyfish'

	input:
	file jellyfish from jellyfiles

	script:
	"""
		jellyfish histo -t 4 ${jellyfish} > ${jellyfish.baseName}.histo
	"""
}

process hifiasm {
	
	conda 'hifiasm'
	memory '128 GB'
	executer 'slurm'

	input:
	file pacbio from files_hifiasm
	tuple val x, tuple file one, file two from files_hic_hifiasm

	script:
	"""
				
	"""
}

workflow {
	 
}