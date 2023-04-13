# variant-ID (v1.0.0)
This Nextflow pipeline can identify variants from next-generation sequencing data for a cohort of samples.

## Contents
- [Workflow](#workflow)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [Example](#example)
- [License](#license)
- [Citations](#citations)

Workflow
--------
The user provides raw paired-end sequencing reads and a reference genome in FASTA format. [**Trim Galore v0.6.10**](https://github.com/FelixKrueger/TrimGalore) is used to perform adapter and quality trimming. [**FastQC v0.11.9**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used to provide some quality control information for both the raw and trimmed reads. [**BBMap v39.01**](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) is used to align reads to the reference genome. Duplicate reads are then marked with [**Picard v2.27.5**](https://broadinstitute.github.io/picard/) and variants are called for each sample using HaplotypeCaller in [**GATK v4.3.0.0**](https://gatk.broadinstitute.org/hc/en-us). Joint genotyping is then performed for the entire cohort and the data are filtered.

Usage
-----
	Please put your raw paired-end FASTQ files into a new folder titled "reads/" and ensure 
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
    
		--ploidy	        Specifies the ploidy of the organism.
					Defaults to "2" for diploid.
						
		--memory		Specifies the amount of RAM available for analysis.
					Defaults to "48G".
						
		--customInterval	Specifies a genomic interval for which variants are called.
					Defaults to the entire genome.
								
		--outputBase		Specifies the basename of the final VCF output file.
					Defaults to "cohort".
							
	Outputs:
		results/    
		|
		|---alignments/		Sequence alignments to the reference genome in Binary
		|			Alignment Map (BAM) format.
		|
		|---fastqc_raw/		Quality control reports for the raw sequencing reads.
		|
		|---fastqc_trimmed/	Quality control reports for the processed sequencing reads.
		|
		|---gvcf/	        Identified variants for individual samples in Genomic 
		|		        Variant Call Format (GVCF).
		|
		|---mapping_stats/	Information regarding the sequence alignments (e.g., coverage).
		|
		|---trimmed_reads/	Processed sequencing reads.
		|
		`---vcf/		The final joint genotyping output in Variant Call Format (VCF).
    
Dependencies 
------------
The following programs must be installed on your computer:
* [**Nextflow**](https://github.com/nextflow-io/nextflow) (v22.10.6 or higher)
* [**Singularity**](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (v3.8.7 or higher)
 
An easy way to install these is with [**Mamba**](https://github.com/mamba-org/mamba). This doesn't require the user to be root so can help to manage software packages on HPC systems.

With Mamba installed, one can fetch the dependencies using the following commands:
````
```
mamba install -c bioconda nextflow
mamba install -c conda-forge singularity
```
````

Containers for software used in the workflow are hosted on [**Docker Hub**](https://hub.docker.com/u/bryoinformatics) and will be downloaded automatically.

Example
------------
A working example can be executed using the sample data distributed with this package.

Step 1: Download the repository and change directory.
````
```
git clone https://github.com/bryoinformatics/variant-ID.git
cd variant-ID
```
````

Step 2a: Make a new directory called "reads/" and move the sample FASTQ files into it. 
Step 2b: Make a new directory called "reference/" and move the sample reference genome into it.
````
```
mkdir reads
cp sample_data/*.fastq.gz reads/
mkdir reference
cp genome.fa.gz* reference/
```
````

Step 3: Run the pipeline.
````
```
nextflow run main.nf
```
````

License
-------
This workflow is released under a GNU General Public License (v3.0).

Citations
---------
* Andrews, S. (2019). FastQC: A quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Broad Institute. (2022). “Picard.” Broad Institute, GitHub repository. http://broadinstitute.github.io/picard/
* Bushnell, B., Rood, J., & Singer, E. (2017). BBMerge – Accurate paired shotgun read merging via overlap. PLOS ONE, 12(10), e0185056. https://doi.org/10.1371/journal.pone.0185056
* Bushnell, B. (2022). "BBMap." SourceForge repository. https://sourceforge.net/projects/bbmap/
* Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008
* Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), Article 4. https://doi.org/10.1038/nbt.3820
* Krueger, F. (2023). Babraham Bioinformatics—Trim Galore! (v0.6.10). https://zenodo.org/badge/latestdoi/62039322
* Kurtzer, G. M., Sochat, V., & Bauer, M. W. (2017). Singularity: Scientific containers for mobility of compute. PLOS ONE, 12(5), e0177459. https://doi.org/10.1371/journal.pone.0177459
* Van der Auwera, G. A., & O’Connor, B. D. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st ed.). O’Reilly Media.

