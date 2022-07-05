#!/bin/bash

# Enabling bash strict mode: (Read more about it at: http://redsymbol.net/articles/unofficial-bash-strict-mode/#expect-nonzero-exit-status)
set -euo pipefail

# Importing configuration
source /MOUNT/config.sh

# uncompressing all fastq.gz files
for i in $(find input/ -name "*fastq.gz" -nowarn )
do echo "uncompressing $i"
pigz -d $i
done
wait

# Making index and output directories, incase there aren't any.
! test -d index && mkdir index
! test -d output && mkdir output

# Building indeces to cache [ Duration: 1.5 hours (on Human genome)]
# ---------------------------------------------------------------------------------------------------------
echo "###################"
echo "##### DASiRe Step 1 of 7: Building indeces for STAR and Kallisto"
echo "###################"
# Building Kallisto's index
if ! test -f /MOUNT/index/kallisto-index
 then echo "Indexing reference genome with Kallisto"
 	kallisto index -i /MOUNT/index/kallisto-index $transcripts_fasta
 else echo "Index found at index/kallisto-index, reusing it:"
 	ls /MOUNT/index/kallisto-index -lah
fi

# Building STAR's index
if ! test -f /MOUNT/index/starindex/genomeParameters.txt
	then mkdir -p /MOUNT/index/starindex || true # to allow mkdir to fail gracefully, we add "|| true"
		echo "Indexing reference genome with STAR"
		STAR \
		--runMode genomeGenerate \
		--genomeDir /MOUNT/index/starindex \
		--genomeFastaFiles $fasta \
		--runThreadN $nCores \
		--sjdbGTFfile $gtf \
		-sjdbOverhang 100 \
		--outFileNamePrefix output/STAR/
	else echo "Indexes found at index/starindex, reusing them:"
	ls /MOUNT/index/starindex/ -lah
 fi 

# Alignments
# ---------------------------------------------------------------------------------------------------------
echo "###################"
echo "##### DASiRe Step 2 of 7: Alignments"
echo "###################"
## Checking if fastq files are paired or single.
for i in $(find input/ -name "*fastq" -nowarn  | sed 's/..fastq$//g' | sed 's/^input\///g'| sort | uniq)
	do fastqcount=$(find input/ -name "${i}?.fastq" | wc -l)
	echo "Checking paires state of sample: input/${i}*  | Number of Fastqfiles: $fastqcount"
	echo "###################"
	
	pairstate="paired"
	if [ $fastqcount -eq 2 ] && [ $pairstate != "single" ] ; then
		if [ -z ${prevcount+set} ] ; then prevcount=${fastqcount}; fi
		if ! [ $prevcount -eq $fastqcount ] ; then echo "DASiRe takes either unpaired or paired reads per session, but not mixed" && exit 1 ; fi
	  echo 'Fastq files found are paired' && pairstate="paired"
	else
  	echo 'Fastq files found are not paired' && pairstate="single"
	fi
done

if [ -z ${adapters+x} ]; then adapters=""; fi
## To use or not to use: GNU Parallel
## A. Parallel grants the capacity to parallelize asyncronous per sample execution. - You cannot be greedy when resources are available.
###or
## B. use $nCores as arguments for each step, allowing parallelisation in a syncronous per task execution. - kills the laptop/PC compute market.
#  We go with A. 
#  Keep the tool functionally within range on laptops and servers alike.

if [ $pairstate == "paired" ]; then
find input/ -name '*fastq' -nowarn  | sed 's/..fastq$//g' | sed 's/^input\///g'| sort | uniq | sed 's/_$//g'| parallel --will-cite -j $nCores "

	echo Transcript quantification with Kallisto for Sample: {}
	! test -d output/pseudocounts/{} && mkdir -p output/pseudocounts/{} || true # to allow mkdir to fail gracefully, we add '|| true'
	if ! test -f output/pseudocounts/{}/abundance.tsv 
		then kallisto quant -i /MOUNT/index/kallisto-index -o output/pseudocounts/{} --bias input/{}_1.fastq input/{}_2.fastq 2> output/pseudocounts/{}/kallisto.stdout
	fi

	echo Genome Alignment with STAR for Sample: {}
	#! test -d output/STAR/ && mkdir -p output/STAR/ || true # to allow mkdir to fail gracefully, we add '|| true'
	if ! test -f output/STAR/{}-Aligned.sortedByCoord.out.bam
		then STAR --genomeDir /MOUNT/index/starindex --readFilesIn input/{}_1.fastq input/{}_2.fastq --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix output/STAR/{}-
	fi
	if ! test -f output/STAR/{}-Aligned.sortedByCoord.out.bam.bai
	then echo indexing {}
		samtools index output/STAR/{}-Aligned.sortedByCoord.out.bam output/STAR/{}-Aligned.sortedByCoord.out.bam.bai
	fi
	"
	fi

if [ $pairstate == "single" ]; then
find input/ -name '*fastq' -nowarn  | sed 's/.fastq$//g' | sed 's/^input\///g'| sort | uniq | parallel --will-cite -j $nCores "echo Sample: {}  Aligning to the genome with STAR \& Kallisto

	echo Transcript quantification with Kallisto for Sample: {}
	! test -d output/pseudocounts/{} && mkdir -p output/pseudocounts/{} || true # to allow mkdir to fail gracefully, we add '|| true'
	if ! test -f output/pseudocounts/{}/abundance.tsv 
		then kallisto quant -i /MOUNT/index/kallisto-index -o output/pseudocounts/{} --bias --single -l $fragment_length -s 20 input/{}.fastq 2> output/pseudocounts/{}/kallisto.stdout
	fi

	echo Genome Alignment with STAR for Sample: {}
	if ! test -f output/STAR/{}-Aligned.sortedByCoord.out.bam
		then STAR --genomeDir /MOUNT/index/starindex --readFilesIn input/{}.fastq --outSAMattributes NH HI AS nM NM MD --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastxm --outFileNamePrefix output/STAR/{}-
	fi
	if ! test -f output/STAR/{}-Aligned.sortedByCoord.out.bam.bai
	then echo indexing {}
		samtools index output/STAR/{}-Aligned.sortedByCoord.out.bam output/STAR/{}-Aligned.sortedByCoord.out.bam.bai
	fi"

fi
wait

# MultiQC for STAR and Kallisto
:> /tmp/multiqc_file_list.txt  # start and clear the file
find /MOUNT/output/pseudocounts/ -name "*.stdout" >> /tmp/multiqc_file_list.txt
find /MOUNT/output/STAR -name "*Log.final.out" >> /tmp/multiqc_file_list.txt
pushd /MOUNT/output/
multiqc -f -l /tmp/multiqc_file_list.txt 
popd
# Analysis with Experiment Design.
# ---------------------------------------------------------------------------------------------------------
echo "###################"
echo "##### DASiRe Step 3 of 7: Majiq AS events"
echo "###################"
		if ! test -f output/MAJIQ/build/splicegraph.sql
			then
		 	mkdir -p /MOUNT/output/MAJIQ/|| true
		
			MajiqConfig=/MOUNT/output/MAJIQ/config.txt
			echo "[info]" > $MajiqConfig
			echo "readlen=75" >> $MajiqConfig
			echo "bamdirs=/MOUNT/output/STAR" >> $MajiqConfig
			echo "genome=hg38" >> $MajiqConfig
			echo "strandness=None" >> $MajiqConfig
			echo "[experiments]" >> $MajiqConfig
			for j in $(let a=$(wc -l /MOUNT/input/metadata.txt | cut -f1 -d ' ')-1 && cat /MOUNT/input/metadata.txt | tail -n $a | cut -f2 | uniq);
			 do echo $j=$(for i in $(grep $j input/metadata.txt | cut -f 1); do basename $(find output/STAR/ -name "$i*bam") .bam ; done)| sed 's/ /,/g'
			done >> $MajiqConfig

			echo "building MAJIQ reference ..."
			majiq build $gff -c $MajiqConfig -j $nCores -o /MOUNT/output/MAJIQ/build
		fi
		let metadata_ln_nos=$(wc -l /MOUNT/input/metadata.txt | cut -f1 -d ' ')-1  # to capture metadata without headers
		if ! test -f output/MAJIQ/voila_results_all.tsv
			then
				for i in $(tail -n $metadata_ln_nos /MOUNT/input/metadata.txt | cut -f2 | uniq);
				do	
					if [ $i != 'control' ]; then
							#-Depreciated command.
							#-get all .majiq files which were created with build
							#-majiqlist=$(ls -1p /MOUNT/output/MAJIQ/build/*.majiq | xargs echo)	
							#-majiq psi $majiqlist -j $nCores -o /MOUNT/output/MAJIQ/$i -n "$i" # $i = 'control\treated'

						#Building up to the majiq deltapsi
						:> /tmp/deltapsi-params.txt
						let k=0 || true #I didn't realize assigning k=0 in bash reads out as 1 in POSIX signal. 
						for j in $(tail -n $metadata_ln_nos /MOUNT/input/metadata.txt | cut -f2 | uniq);
							do let k=$k+1; echo "-grp${k} $(for l in $(grep $j input/metadata.txt | cut -f 1); do find output/MAJIQ -name "$l*majiq" && echo " " ; done) "
						done >> /tmp/deltapsi-params.txt
						echo ' -n' >> /tmp/deltapsi-params.txt
						for j in $(tail -n $metadata_ln_nos /MOUNT/input/metadata.txt  | cut -f2 | uniq);
			 				do echo " $j" 
						done >> /tmp/deltapsi-params.txt
						majiq deltapsi -j $nCores -o /MOUNT/output/MAJIQ/deltapsi $(cat /tmp/deltapsi-params.txt | tr -d '\n')

					# create voila.tsv outputfiles
					voila tsv /MOUNT/output/MAJIQ/build/splicegraph.sql /MOUNT/output/MAJIQ/deltapsi/*.voila -f /MOUNT/output/MAJIQ/voila_results.tsv 
					voila tsv /MOUNT/output/MAJIQ/build/splicegraph.sql /MOUNT/output/MAJIQ/deltapsi/*.voila -f /MOUNT/output/MAJIQ/voila_results_all.tsv --show-all
					break
					fi					
				done
			else echo "MAJIQ looks like it ran. Please take a look at the output files and remove them, if it is not the run you wanted. And run the pipeline again."
		fi
echo "###################"
echo "##### DASiRe is running in R from here."
echo "###################"

Rscript /deseq2_dexseq_isoformswitchanalyzer.Rscript -b /MOUNT/output/STAR/ -l $pairstate --gtffile $gtf -o /MOUNT/output/ -m $metadata -f $transcripts_fasta -p $nCores