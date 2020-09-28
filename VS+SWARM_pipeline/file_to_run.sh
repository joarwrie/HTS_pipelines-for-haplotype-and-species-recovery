#!/bin/bash

# Require VSEARCH-2.14.1
# Require cutadapt-2.8
# Require swarm-3.0.0
# Require usearch-11.0.667 (32 bits version)
# Require obitools-1.2.11

# Demultiplexing and primer removal
for i in $(ls *R1_001.fastq.gz);do
	Nom=$(echo $i | cut -f1 -d"_")
	reverse=$(basename $i "R1_001.fastq.gz")"R2_001.fastq.gz"
	cutadapt -g file:barcode_forward.fasta -G file:barcode_reverse.fasta -e 0.15 --no-indels -o ${Nom}_{name1}_{name2}_R1_cut.fastq -p ${Nom}_{name1}_{name2}_R2_cut.fastq $i $reverse
done

# Merging read pairs
for i in $(ls *ext*PCR*R1_cut.fastq);do
	Nom=$(echo $i | cut -f1,2,3 -d"_")
	reverse=$(basename $i R1_cut.fastq)"R2_cut.fastq"
	vsearch --fastq_mergepairs "$i" --reverse $reverse --fastqout ${Nom}_merged.fastq --relabel_keep --relabel "${Nom};read" --fastq_maxdiffs 100 --fastq_maxee 2 --fastq_minmergelen 402 --fastq_maxmergelen 422 --eeout
done

# Concatenation of all files and dereplication (removing singletons)
cat *_merged.fastq > AllSamples_merged_final.fastq
vsearch --derep_fulllength AllSamples_merged_final.fastq --output AllSamples_dereplic.fasta --sizeout --minuniquesize 2 --notrunclabels
vsearch --fastq_filter AllSamples_merged_final.fastq --fastaout AllSamples_merged_final.fasta

# Remove sequences with Ns 
obigrep --fasta --nuc --uppercase -s '^[ATCG]+$' AllSamples_dereplic.fasta > AllSamples_dereplic_filtered.fasta
cat AllSamples_dereplic_filtered.fasta | sed 's/>[^ ]*   \(.*\)/>\1/g' > AllSamples_dereplic_filtered_modif.fasta

# Clustering
swarm -d1 -f -b200 -w AllSamples_centroids.fasta -o AllSamples_OTUlist.txt -z AllSamples_dereplic_filtered_modif.fasta

# Creating contingency table
usearch -otutab AllSamples_dereplic_withinSamples.fasta -otus AllSamples_dereplic_filtered_modif.fasta -otutabout AllSamples_finalTAB.csv -id 1 -notmatched AllSamples_unmapped.fasta

# Filtering for index-jumping and replicates
Rscript Data_treatment.R

