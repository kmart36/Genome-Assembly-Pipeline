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
        path "${sample.baseName}.jf"

	script:
	"""
		jellyfish count -C -m 21 -s 16000000000 -t 10 ${sample} -o ${sample.baseName}.jf

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
	executer 'slurm'
	cpus 4

	input:
	file pacbio
	file hic1
	file hic2

	output:
	path "${pacbio.baseName}.asm.hic.p_ctg.gfa"	

	script:
	"""
	hifiasm -l 0 -o ${pacbio.baseName}.asm -t 32 --h1 ${hic1} --h2 ${hic2} ${pacbio}	
	"""
}

process gfa2fasta {
	
	executer 'slurm'
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
		 
	input:
	file fasta from asm_fasta	
	
	script:
	"""
	bwa index -p ${fasta.baseName}_bwa_index ${fasta}
	"""
}

process bwa_mem {
	
	input:
	file val

	script:
	"""
	bwa mem -t 4 _bwa_index /home/kam071/SampleReads/Rhagonycha_fulva/ERR6054569_1.fastq /home/kam071/SampleReads/Rhagonycha_fulva/ERR6054569_2.fastq 2> bwa.err > bwa_fulva_HiC.sam
	"""
}

workflow {
	 pacBio = Channel.fromPath('/home/kam071/SampleReads/Photinus_scintillans/RawHiFiData/**/*.fastq*')
	 hic_1 = Channel.fromPath('/home/kam071/SampleReads/Photinus_scintillans/RawHiCData/*_R1.fastq*')
	 hic_2 = Channel.fromPath('/home/kam071/SampleReads/Photinus_scintillans/RawHiCData/*_R2.fastq*')
	 fastQC(pacBio)
	 jellyfish(pacBio)
	 genomescope(jellyfish.out)
	 hifiasm(pacBio, hic_1, hic_2)
	 gfa2fasta(hifiasm.out)
}