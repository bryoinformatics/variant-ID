#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*---------------------------------------------------------------------------------------------

	Please put your raw paired-end FASTQ files into the folder "reads/" and ensure the 
		that the files end with the suffix ".fastq.gz". The sample name should precede
		the read pair  designation. e.g., the files for SampleA would be 
		"SampleA_1.fastq.gz" and "SampleA_2.fastq.gz".
		
	Please put your reference genome in FASTA format in the folder "reference/". This file
		should be compressed and indexed with bgzip. Ensure that the file name is 
		"genome.fa.gz".
		
	Usage:
		nextflow run main.nf
	
	Inputs:
	
		(Optional)
		
		--ploidy			Specifies the ploidy of the organism.
								Defaults to "2" for diploid.
						
		--memory			Specifies the amount of RAM available for analysis.
								Defaults to "48G".
						
		--customInterval	Specifies a genomic interval for which variants are called.
								Defaults to the entire genome.
								
		--outputBase		Specifies the basename of the final VCF output file.
								Defaults to "cohort".
							
							
	Outputs:
		results/    
		|
		|---alignments/			Sequence alignments to the reference genome in Binary
		|							Alignment Map (BAM) format.
		|
		|---fastqc_raw/			Quality control reports for the raw sequencing reads.
		|
		|---fastqc_trimmed/		Quality control reports for the processed sequencing reads.
		|
		|---gvcf/				Identified variants for individual samples in Genomic 
		|							Variant Call Format (GVCF).
		|
		|---mapping_stats/		Information regarding the sequence alignments (e.g., coverage).
		|
		|---trimmed_reads/		Processed sequencing reads.
		|
		`---vcf/				The final joint genotyping output in Variant Call Format (VCF).
		
---------------------------------------------------------------------------------------------*/		
 
 
/*---------------------------
	Define the parameters
---------------------------*/


params.genome = "./reference/genome.fa.gz"
params.reads = "reads/*_{1,2}.fastq.gz"
params.ploidy = 2
params.memory = "48G"
params.customInterval = false
params.outputBase = "cohort"


/*-------------------------
	Define the workflow
-------------------------*/


workflow {
	reads_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )
	FASTQC( reads_ch )
	trim_ch = TRIM( reads_ch )
	FASTQC_TRIM( trim_ch )
	genome_ch = Channel.fromPath( params.genome, checkIfExists:true )
	INDEX( genome_ch )
	ALIGN( trim_ch, genome_ch )
	CALL( ALIGN.out.bam_files, genome_ch ) | collect | CAT
	GENOTYPE( genome_ch, CAT.out.cohort, CAT.out.interval)
	
}


/*-------------------------------------------------
	Define the processes called in the workflow
-------------------------------------------------*/


process FASTQC {
	tag "Processing $sampleId"
	
    publishDir './results/fastqc_raw', mode:'copy'

	input:
		tuple val(sampleId), file(reads_ch)
	
	output:
		tuple val(sampleId), path("${sampleId}_{1,2}_fastqc.html")
		tuple val(sampleId), path("${sampleId}_{1,2}_fastqc.zip")
	
    script:
		"""
		fastqc ${reads_ch} -q -t $task.cpus
		"""
}


process TRIM {
	tag "Processing $sampleId"
	
	publishDir './results/trimmed_reads/', pattern: '*_val_*.fq.gz', mode:'copy'

	input:
		tuple val(sampleId), file(reads_ch)
	
	output:
		tuple val(sampleId), path("${sampleId}_{1,2}_val_{1,2}.fq.gz"), emit: trimmed_reads
	
    script:
		"""
		trim_galore --max_n 4 --cores $task.cpus --paired --gzip ${reads_ch}
		"""
}


process FASTQC_TRIM {
	tag "Processing $sampleId"
	
    publishDir './results/fastqc_trimmed', mode:'copy'

	input:
		tuple val(sampleId), path(trimmed_reads)
	
	output:
		tuple val(sampleId), path("${sampleId}_{1,2}_val_{1,2}_fastqc.html")
		tuple val(sampleId), path("${sampleId}_{1,2}_val_{1,2}_fastqc.zip")
	
    script:
		"""
		fastqc ${trimmed_reads[0]} ${trimmed_reads[1]} -q -t $task.cpus
		"""
}


process INDEX {
	tag "Indexing reference genome"
	
	publishDir './reference', mode:'copy'
	
	input:
		path genome_ch
		
	output:
		path "${genome_ch}.fai"
		path "${genome_ch.SimpleName}.dict"
		
    script:
		"""
		samtools faidx ${genome_ch}
		gatk CreateSequenceDictionary \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-O ${genome_ch.SimpleName}.dict
		"""
}
	
	
process ALIGN {
	tag "Processing $sampleId"
	
    publishDir './results/mapping_stats', mode:'copy', pattern:'*.{coverage,log}'
	publishDir './results/alignments', mode:'copy', pattern:'*.bam'

	input:
		tuple val(sampleId), path(trimmed_reads)
		each genome_ch
	
	output:
		path("${sampleId}.coverage")
		path("${sampleId}.log")
		tuple val(sampleId), path("${sampleId}.bam"), emit: bam_files
	
    script:
		"""
		bbmap.sh -Xmx${params.memory} ref=${genome_ch} in=${trimmed_reads[0]} \
			in2=${trimmed_reads[1]} out=${sampleId}.sam covstats=${sampleId}.coverage \
			bamscript=${sampleId}_bs.sh |& tee ${sampleId}.log
		sh ${sampleId}_bs.sh
		picard MarkDuplicates -I ${sampleId}_sorted.bam -O ${sampleId}_duplicates.bam \
				-M ${sampleId}_duplicates.metrics
		samtools addreplacerg -r "ID:foo" -r "SM:${sampleId}"  -o ${sampleId}_gatk.bam \
			${sampleId}_duplicates.bam
		samtools sort -o ${sampleId}.bam ${sampleId}_gatk.bam
		"""
}


process CALL {
	tag "Processing $sampleId"
	
	publishDir './results/gvcf', mode:'copy'
	
	input:
		tuple val(sampleId), path(bam_files)
		each genome_ch
	
	output:
		path("${sampleId}.gvcf*"), emit: gvcf_files
		
    script:
	if ( params.customInterval == false )
		"""
		samtools index ${sampleId}.bam
		gatk HaplotypeCaller --java-options "-Xmx${params.memory} -Xms${params.memory}" \
			-R ${genome_ch} -I ${sampleId}.bam -O ${sampleId}.gvcf -ERC GVCF \
			-ploidy ${params.ploidy}
		"""
	else
		"""
		samtools index ${sampleId}.bam
		gatk HaplotypeCaller --java-options "-Xmx${params.memory} -Xms${params.memory}" \
			-R ${genome_ch} -I ${sampleId}.bam -O ${sampleId}.gvcf -ERC GVCF \
			-ploidy ${params.ploidy} -L ${params.customInterval}
		"""
}


process CAT {
	tag "Processing $sampleId"
	
	input:
		path(gvcf_files)
		
	output:
		path("cohort.txt"), emit: cohort
		path("intervals.bed"), emit: interval
		
    shell:
		$/
		ls "$PWD/"results/gvcf/*.gvcf > cohort_temp1.txt
		cat cohort_temp1.txt | sed 's/\.gvcf//g' > cohort_temp2.txt
		sed -i 's/.*results\/gvcf\///g' cohort_temp2.txt
		paste cohort_temp2.txt cohort_temp1.txt >cohort.txt
		awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' \
			"$PWD/"reference/*.fai > intervals.bed
		/$
}


process GENOTYPE {
	tag "Genotyping cohort"
	
    publishDir './results/vcf', mode:'copy'

	input:
		path genome_ch
		path cohort
		path interval
		
	output:
		path("${params.outputBase}.vcf.gz*")
		
    script:
	if ( params.customInterval == false )
		"""
		samtools faidx ${genome_ch}
		gatk CreateSequenceDictionary \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-O ${genome_ch.SimpleName}.dict
		gatk GenomicsDBImport \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			--sample-name-map ${cohort} --genomicsdb-workspace-path database -L ${interval}
		gatk GenotypeGVCFs \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-V gendb://database -O gatk.vcf -L ${interval}
		gatk VariantFiltration \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-V gatk.vcf -O gatk.filtered.vcf \
			--filter-name "QD_filter" --filter-expression "QD < 2.0" \
			--filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
			--filter-name "FS_filter" --filter-expression "FS > 60.0" \
			--filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
			--filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
			--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
			-L ${interval}
		bgzip -c gatk.filtered.vcf >gatk.filtered.vcf.gz
		bcftools index gatk.filtered.vcf.gz
		bcftools view -f 'PASS' gatk.filtered.vcf.gz -O v -o gatk.pass.vcf
		gatk SelectVariants \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-V gatk.pass.vcf --select-type-to-include SNP --select-type-to-include INDEL \
			-O gatk.final.vcf
		bgzip -c gatk.final.vcf >${params.outputBase}.vcf.gz
		bcftools index ${params.outputBase}.vcf.gz
		"""
	else
		"""
		samtools faidx ${genome_ch}
		gatk CreateSequenceDictionary \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-O ${genome_ch.SimpleName}.dict
		gatk GenomicsDBImport \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			--sample-name-map ${cohort} --genomicsdb-workspace-path database \
			-L ${params.customInterval}
		gatk GenotypeGVCFs \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-V gendb://database -O gatk.vcf -L ${params.customInterval}
		gatk VariantFiltration \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-V gatk.vcf -O gatk.filtered.vcf \
			--filter-name "QD_filter" --filter-expression "QD < 2.0" \
			--filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
			--filter-name "FS_filter" --filter-expression "FS > 60.0" \
			--filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
			--filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
			--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
			-L ${params.customInterval}
		bgzip -c gatk.filtered.vcf >gatk.filtered.vcf.gz
		bcftools index gatk.filtered.vcf.gz
		bcftools view -f 'PASS' gatk.filtered.vcf.gz -O v -o gatk.pass.vcf
		gatk SelectVariants \
			--java-options "-Xmx${params.memory} -Xms${params.memory}" -R ${genome_ch} \
			-V gatk.pass.vcf --select-type-to-include SNP --select-type-to-include INDEL \
			-O gatk.final.vcf
		bgzip -c gatk.final.vcf >${params.outputBase}.vcf.gz
		bcftools index ${params.outputBase}.vcf.gz
		"""
}