#!/bin/bash

# Require OBITOOLS-1.2.11
# Require swarm-3.0.0
# Require usearch11.0.667 (32 bits version)

# Merging read pairs
for i in $(ls *R1_001.fastq);do
	Nom=$(echo $i | cut -f1 -d"_")
	reverse=$(basename $i "R1_001.fastq")"R2_001.fastq"
	illuminapairedend $i -r $reverse --sanger --fastq-output --uppercase --score-min=40 | obigrep -p 'mode!="joined"' --sanger --nuc --uppercase > ${Nom}_paired.fastq
done

# Primer removal and demultiplexing on tags to separate extraction and PCR replicates
for i in $(ls *_paired.fastq);do
	Nom=$(echo $i | cut -f1 -d"_")
	ngsfilter -e 2 --sanger --nuc --fastq-output --uppercase -t barcode.txt -u ${Nom}_unidentified.fastq $i > ${Nom}_demultiplexed.fastq
done

# Creation of a unique file for the whole dataset
touch AllSamples_demultiplexed.fastq
for i in $(ls *_demultiplexed.fastq);do
	Nom=$(echo $i | cut -f1 -d"_")
	obiannotate --sanger --fastq-output --uppercase -S Pot:${Nom} $i | obiannotate --sanger --fastq-output --uppercase -S 'Echantillon:sequence["Pot"]+"_"+sequence["sample"]' >> AllSamples_demultiplexed.fastq
done

# Dereplication
obiuniq -m Echantillon --sanger --nuc AllSamples_demultiplexed.fastq > AllSamples_dereplic.fasta

# Length selection and singletons removal
obigrep -l 402 -L 422 --fasta --nuc --uppercase AllSamples_dereplic.fasta | obigrep -p 'sequence["count"]!=1' --fasta --nuc --uppercase > AllSamples_length.fasta

# Formatting file headers
cat AllSamples_length.fasta | sed 's/>\([^ ]*\) .*count=\([0-9]*\); .*/>\1;size=\2/g' > AllSamples_formatted.fasta

# Clustering
swarm -d1 -f -b200 -w AllSamples_centroids.fasta -o AllSamples_OTUlist.txt -z AllSamples_formatted.fasta

# Creation of the distribution table
Tableau_info_seq.py AllSamples_length.fasta Tab_final_distri.txt

# Identification of the maximum number of reads in a index control sample and filtering for index control and replicates
Rscript Data_treatment.R
