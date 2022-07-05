#!/bin/bash

# Filepaths that begin with '/' suggests that these are absolute filepaths, which are needed for processes within the docker.
# /MOUNT within the docker  == DASiRe/MOUNT/ directory in your local machine.

#This file hosts the variables for DASiRe's preprocessing.
nCores=4
fasta=/MOUNT/input/GRCh38.primary_assembly.genome.fa
gtf=/MOUNT/input/gencode.v36.annotation.gtf
gff=/MOUNT/input/gencode.v36.annotation.gff3
transcripts_fasta=/MOUNT/input/gencode.v36.transcripts.fa
metadata=/MOUNT/input/metadata.txt


#optional
fragment_length=170
