#!/bin/bash

# Demultiplexing on the set of tags used for PCR and extraction replicates identification
# Primer removal

# Require cutadapt-2.3

for i in $(ls *R1_001.fastq.gz);do
	Nom=$(echo $i | cut -f1 -d"_")
	reverse=$(basename $i "R1_001.fastq.gz")"R2_001.fastq.gz"
	cutadapt -g file:barcode_forward.fasta -G file:barcode_reverse.fasta -e 0.15 --no-indels -o ${Nom}_{name1}_{name2}_R1_cut.fastq -p ${Nom}_{name1}_{name2}_R2_cut.fastq $i $reverse
done

mkdir Demult_files
mv *ext*_PCR*.fastq Demult_files/

# Running DADA2

Rscript DADA2.R
