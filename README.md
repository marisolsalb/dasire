# DASiRe

## What is DASiRe?

**Direct Alternative Splicing Regulator predictor (DASiRe)** is a web application that allows non-expert users to perform different types of splicing analysis from RNA-seq experiments and also incorporates ChIP-seq data of a DNA-binding protein of interest to evaluate whether its presence is associated with the splicing changes detected in the RNA-seq dataset. 

DASiRe is an accessible web-based platform that performs the analysis of raw RNA-seq and ChIP-seq data to study the relationship between DNA-binding proteins and alternative splicing regulation. It provides a fully integrated pipeline that takes raw reads from RNA-seq and performs extensive splicing analysis by incorporating the three current methodological approaches to study alternative splicing: isoform switching, exon and event-level. Once the initial splicing analysis is finished, DASiRe performs ChIP-seq peak enrichment in the spliced genes detected by each one of the three approaches. 

## Serverside
Start the serverside app by running the following command from the top directory of this git.
(server-update.sh is a prerequisite for internal use only)

```
 docker run -v $(pwd)/serverside/R:/srv/shiny-server/ --rm --name 'dasire' -p 3838:3838 -d marisalb/dasire:server
```

## Preprocessing 
Start preprocessing pipeline by first dropping all of your inputs in MOUNT/input directory and re-write your config file.

```
MOUNT/
├──input/
|   ├── adapters_pe.fa
|   ├── ENCFF239POR_1.fastq
|   ├── ENCFF239POR_2.fastq
|   ├── ENCFF250AJS_1.fastq
|   ├── ENCFF250AJS_2.fastq
|   ├── ENCFF301JRH_1.fastq
|   ├── ENCFF301JRH_2.fastq
|   ├── ENCFF803NZJ_1.fastq
|   ├── ENCFF803NZJ_2.fastq                     # Also, mention these files in the config file
|   ├── gencode.v36.annotation.gff3             # gff=/MOUNT/input/gencode.v36.annotation.gff3
|   ├── gencode.v36.annotation.gtf              # gtf=/MOUNT/input/gencode.v36.annotation.gtf
|   ├── gencode.v36.transcripts.fa              # transcripts_fasta=/MOUNT/input/gencode.v36.transcripts.fa
|   ├── GRCh38.primary_assembly.genome.fa       # fasta=/MOUNT/input/GRCh38.primary_assembly.genome.fa
|   └── metadata.txt                            # metadata=/MOUNT/input/metadata.txt
└── config.sh                          <---# Hosted above the input directory.
```
Also add the number of cores you want DASiRe to use with `nCores=4`. I'd recommend 4 on a laptop with 6-8 cores.

The pipeline uses the following software versions:
- GNU parallel 20161222
- MultiQC v1.12
- Kallisto:0.48.0
- STAR: 2.7.10a
- samtools Version: 1.7 (using htslib 1.7-2)
- majiq 2.1
- genomicranges:1.46.1
- genomicfeatures:1.46.1
- genomicalignments:1.30.0
- isoformswitchanalyzer:1.16.0
- deseq2:1.34.0
- dexseq:1.40.0
- rsubread:2.8.1

To start the DASiRe pipeline, please run the following command in this git's main directory.
`docker run -v $(pwd)/MOUNT:/MOUNT --user $(id -u):$(id -g) --name 'dasire-preprocessing' --rm marisalb/dasire:client`