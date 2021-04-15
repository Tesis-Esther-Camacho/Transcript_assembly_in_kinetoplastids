# Transcript Asembly in Kinetoplastids


## INTRODUCTION

This pipeline can be used in order to assembly the transcriptomes of kinetoplastids from RNA-sequencing data. It has been successfully tested on different Leishmania species: *L. infantum*, *L. major*, *L. donovani* and *L. braziliensis*.

## REQUIREMENTS

This pipeline requires the following dependences:

* [Python 2.7.X](https://www.python.org/downloads/)
* [Bowtie-2.2.9](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools-1.2](http://samtools.sourceforge.net/)
* [Cufflinks-2.0.2](http://cole-trapnell-lab.github.io/cufflinks/) 
* [Blastn 2.2.28+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

All these softwares must be in the local PATH.

## INSTALATION

Make sure you have installed all the dependencies before try to run the pipeline!

Main script (found in the downloads folder):

* kineto_v5.1.py

The downloaded folder contains all the in-house additional scripts needed to do the complete pipeline:

* leader_pair_v1.py
* leader_single_v1.py
* polyA_pair_v1.py
* polyA_single_v1.py
* SL_cut_v4.py
* SL_recut_v1.py
* polyA_cut_v4.py
* polyA_recut_v1.py
* intergenic_SL.py
* intergenic_polyA.py
* SL_assign.py
* multi_read_polyA_assign.py
* polyA_re-assign.py
* polycistronic_transcripts_v2.py
* truncated_transcripts_v2.py

## USAGE

```	
kineto.py <step> [options] <organism_name>

Input files:
	
	organism_name.fastq (phred 33 by default)
		organism_name_1.fastq, organism_name_2.fastq for paired-end reads (--pair)
	organism_name.fasta
	organism_name.gff (from triTrypDB)
		
Transcriptome assembling steps (organism_v0):
### STEP 1:
	-1	SeqLider search
		Polya search 
		Initial alignment (Bowtie2)
		Primary assembly (Cufflinks)
		SL division 

Before proceeding with STEP 2 a manual review is needed. This manual review consists on asigning the strand of the transcripts and divide them into smaller transcript if necessary. The file to review is  transcripts_polyA_both_reviewed.gtf

### STEP 2:
	
	-2	PolyA division

Before proceeding with STEP 2 a manual review is needed. This manual review consists on asigning the strand of the transcripts and divide them into smaller transcript if necessary. The file to review is transcripts_SL_both_reviewed.gtf.

### STEP 3:
	
	-3	Initial annotation (Blastn)
		Refinement
		Re-annotation (Blastn)
		Quantification (Cufflinks)

	Note: To execute advanced steps, previous ones must be already done

Options (for transcriptome assembly, no for corrections):

	--sp	Kinetoplastid specie (L=leishmania (default), T=Trypanosoma)
	--pair	For paired-end reads (default single end reads are expected)
	--quals	phred33/phred64 (default = phred33)
	--nt1	mininum nucleotides for SL identifications (default = 8)
	--nt2	mininum nucleotides for polyA identifications (default = 6)

Corrections (put input file in gtf format in transcriptome folder):

	The correction steps are incompatible with transcriptome steps.
	--correct1=<input_file.gtf>	quantification
	--correct2=<input_file.gtf>	annotation + quantification

	Note: The transcriptome must be already assembled in a previous run of the program, 
	and in both cases a new version will be created
```

## AUTHORS

* Alberto Rastrojo Lastras (CBMSO) - arastrojo@cbm.csic.es

## MAINTAINERS

* Alberto Rastrojo Lastras (CBMSO) - arastrojo@cbm.csic.es
* Esther Camacho Cano (CBMSO) - ecamacho@cbm.csic.es
