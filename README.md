# superQ

superQ is a pipeline built for analysis  of PEM-seq data. PEM-seq is a next-generation sequencing method to quantify editing efficiency of CRISPR systems by introducing a random molecular barcode systems.

## Getting Started

superQ consists of four parts:

1.Reads Preprocessing

2.a pipeline to call breakpoint of translocation (https://github.com/robinmeyers/transloc_pipeline)

3.Deduplication according to random molecular barcode

4.find small indels around cut site

### Prerequisites

Please be sure to install: 

1. cutadapt (http://cutadapt.readthedocs.io/en/stable/)

2. fastq-multx (https://github.com/brwnj/fastq-multx)

3. Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

4. macs2 (https://pypi.org/project/MACS2/)

### Installing

```
git clone https://github.com/liumz93/superQ 
```


## Running superQ

### Samples demltiplex
You can use FastMultx.py to demultiplex samples from a library

```
FastMultx.py -i INDEX -1 FASTQ_R1 -2 FASTQ_R2
```

Index file must be look like:

```
sample1	AGCG
sample2	GCCT
sample3	AGGA
sample4	TCAG
```
### Running superQ

```
superQ.py -i inputdir -m meta.txt
```
meta.txt file must be look like:

```
Library	Assembly	Chr	Start	End	Strand	MID	Primer	Adapter	Description
sample1	hg38	chr11	36573328	36573417	-	AGCG	AGGATCTCACCCGGAACAGC	CCACGCGTGCTCTACA	RAG1A_as_bait_site
sample2	hg38	chr11	36573328	36573417	-	GCCT	AGGATCTCACCCGGAACAGC	CCACGCGTGCTCTACA	RAG1A_as_bait_site
sample3	hg38	chr11	36573328	36573417	-	AGGA	AGGATCTCACCCGGAACAGC	CCACGCGTGCTCTACA	RAG1A_as_bait_site
sample4	hg38	chr11	36573328	36573417	-	TCAG	AGGATCTCACCCGGAACAGC	CCACGCGTGCTCTACA	RAG1A_as_bait_site
```


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

This pipeline is developed based on transloc\_pipeline(https://github.com/robinmeyers/transloc_pipeline).