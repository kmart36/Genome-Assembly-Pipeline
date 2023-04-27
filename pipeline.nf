#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workDir = '/home/kam071/SampleReads/Photinus_scintillans'

process fastQC {
	conda 'fastqc'
	executor 'slurm'
	memory '8 GB'

	publishDir '/SampleReads/Photinus_scintillans'	

	input:
	file sample
	
	script:
	"""
		mkdir -p ${sample.baseName}_fastqc
		fastqc --outdir ${sample.baseName}_fastqc ${sample}
	"""
}

process jellyfish {
	
	conda 'jellyfish'
	cpus 8
	memory '100 GB'
	executor 'slurm'

        input:
	file sample
	
	output:
        path "${sample.baseName}_31_IND2.jf"

	script:
	"""
		jellyfish count -C -m 31 -s 16000000000 -t 10 ${sample} -o ${sample.baseName}_31_IND2.jf

	"""
}

process genomescope {
	
	conda 'jellyfish'
	executor 'slurm'
	cpus 1
	memory '4 GB'

	input:
	file jellyfish

	script:
	"""
		jellyfish histo -t 4 ${jellyfish} > ${jellyfish.baseName}.histo
	"""
}

process hifiasm {
	
	conda 'hifiasm'
	memory '128 GB'
	executor 'slurm'
	cpus 4

	input:
	file pacbio
	tuple val(sample_id), path(hic)

	output:
	path "${pacbio.baseName}.asm.hic.p_ctg.gfa"	

	script:
	"""
	hifiasm -l 0 -o ${pacbio.baseName}.asm -t 32 --h1 ${hic[0]} --h2 ${hic[1]} ${pacbio}	
	"""
}

process gfa2fasta {
	
	executor 'slurm'
	cpus 1
	memory '16 GB'	

	input:
	file graph

	output:
	path "${graph.baseName}.fa"

	script:
	"""
	gfatools gfa2fa ${graph} > ${graph.baseName}.fa
	"""
}

process bwa_index {
		 
	executor 'slurm'
	cpus 2
	memory '8 GB'

	input:
	file fasta
	
	output:
	path '*_bwa_index'	
	
	script:
	"""
	bwa index -p ${fasta.baseName}_bwa_index ${fasta}
	"""
}

process bwa_mem {
	
	input:
	file index
	tuple val(sample_id), path(hic)

	script:
	"""
	bwa mem -t 4 _bwa_index ${hic[0]} ${hic[1]} 2> bwa.err > ${index.baseName}.sam
	"""
}

workflow pipeline {
	 pacBio = Channel.fromPath('/home/kam071/SampleReads/Photinus_scintillans/RawHiFiData/**/*.fastq')
	 hic = Channel.fromFilePairs('/home/kam071/SampleReads/Photinus_scintillans/RawHiCData/*_R{1,2}.fastq*', checkIfExists: true)
	 fastQC(pacBio)
         jellyfish(pacBio)
         genomescope(jellyfish.out)
         hifiasm(pacBio, hic)
         gfa2fasta(hifiasm.out)
         bwa_index(gfa2fasta.out)
         bwa_mem(bwa_index.out, hic)
}

workflow simple {
	 pacBio = Channel.fromPath('/home/kam071/SampleReads/Photinus_scintillans/RawHiFiData/**/*.fasta')
	 jellyfish(pacBio)
         genomescope(jellyfish.out)
}

workflow {
	 simple()	 
}