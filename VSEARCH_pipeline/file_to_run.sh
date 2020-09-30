#!/bin/bash

# Require VSEARCH-2.14.1
# Require cutadapt-2.8

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
	vsearch --fastq_mergepairs $i --reverse $reverse --fastqout ${Nom}_merged.fastq --relabel_keep --relabel ${Nom}_read --fastq_maxdiffs 100 --fastq_maxee 2 --fastq_minmergelen 402 --fastq_maxmergelen 422 --eeout
done

# Concatenation of all files and dereplication (removing singletons)
cat *_merged.fastq > AllSamples_merged_final.fastq
vsearch --derep_fulllength AllSamples_merged_final.fastq --output AllSamples_dereplic.fasta --sizeout --minuniquesize 2 --notrunclabels
vsearch --fastq_filter AllSamples_merged_final.fastq --fastaout AllSamples_merged_final.fasta

# Clustering
vsearch --cluster_smallmem AllSamples_dereplic.fasta --id 0.995 --iddef 2 --sizein --sizeout --uc AllSamples_clusteringTAB.csv --centroids AllSamples_centroids.fasta --clusterout_id --usersort

# Mapping reads to OTUs
vsearch -usearch_global AllSamples_merged_final.fasta -db AllSamples_centroids.fasta -otutabout AllSamples_finalTAB.csv -id 0.995 -iddef 2 -notmatched AllSamples_unmapped.fasta

# Filtering for index-jumping and replicates
Rscript Data_treatment.R

