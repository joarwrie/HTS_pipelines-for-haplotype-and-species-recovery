# Require mothur-1.42.0
# Require vsearch-2.9.1

# Demultiplexing and merging pairs of reads
mothur "#make.file(type=fastq, numcols=3, prefix=AllSamples, inputdir=Raw_data/, outputdir=Raw_data/)"
mothur "#make.contigs(file=AllSamples.files, inputdir=Raw_data/, format=illumina1.8+, oligos=Barcode_file_mothur.txt, bdiffs=1, pdiffs=2, processors=16, outputdir=., maxee=2)"
mothur "#screen.seqs(inputdir=., outputdir=., fasta=AllSamples.trim.contigs.fasta, group=AllSamples.contigs.groups, contigsreport=AllSamples.contigs.report, minlength=402, maxlength=422, maxn=0, minoverlap=126, mismatches=10, processors=16)"

# Dereplication
mothur "#unique.seqs(fasta=AllSamples.trim.contigs.good.fasta, format=name)"
mothur "#count.seqs(name=AllSamples.trim.contigs.good.names, group=AllSamples.contigs.good.groups)"

# All reads must be aligned for further steps so we aligned them against the reference database used for species assignment
mothur "#align.seqs(fasta=AllSamples.trim.contigs.good.unique.fasta, reference=RefDB_ascidians_align.fasta)"
mothur "#summary.seqs(fasta=AllSamples.trim.contigs.good.unique.align)"
mothur "#screen.seqs(fasta=AllSamples.trim.contigs.good.unique.align, count=AllSamples.trim.contigs.good.count_table, summary=AllSamples.trim.contigs.good.unique.summary, start=22, end=439, maxhomop=8)"

# Chimera removal using the vsearch tool
mothur "#chimera.vsearch(vsearch=/usr/local/genome2/vsearch-2.9.1/bin/vsearch, fasta=AllSamples.trim.contigs.good.unique.good.align, count=AllSamples.trim.contigs.good.good.count_table, dereplicate=f)"
mothur "#remove.seqs(count=AllSamples.trim.contigs.good.good.count_table, accnos=AllSamples.trim.contigs.good.unique.good.denovo.vsearch.accnos)"
mothur "#remove.seqs(fasta=AllSamples.trim.contigs.good.unique.good.align, accnos=AllSamples.trim.contigs.good.unique.good.denovo.vsearch.accnos)"

# Remove singletons
cat AllSamples.trim.contigs.good.good.pick.count_table | cut -f1,2 | tr "\t" "%" | egrep "%1$" | sed 's/%1//g' > AllSamples.singletons.accnos
mothur "#remove.seqs(count=AllSamples.trim.contigs.good.good.pick.count_table, accnos=AllSamples.singletons.accnos)"
mothur "#remove.seqs(fasta=AllSamples.trim.contigs.good.unique.good.pick.align, accnos=AllSamples.singletons.accnos)"

# obidistribute + muscle pour aligner les séquences au format fasta
mothur "#dist.seqs(fasta=AllSamples.trim.contigs.good.unique.good.pick.pick.align, countends=F, cutoff=0.03)"
mothur "#cluster(column=AllSamples.trim.contigs.good.unique.good.pick.pick.dist, count=AllSamples.trim.contigs.good.good.pick.pick.count_table, method=opti, cutoff=0.005)"

# Pipeline maison pour obtenir le tableau de contingence et le fasta des séquences de référence pour chaque OTU
MakeShared_custom.sh AllSamples.trim.contigs.good.good.pick.pick.count_table AllSamples.trim.contigs.good.names
mothur "#get.oturep(column=AllSamples.trim.contigs.good.unique.good.pick.pick.dist, list=AllSamples.trim.contigs.good.unique.good.pick.pick.opti_mcc.list, name=AllSamples.trim.contigs.good.names)"
Tableau_info_seq.py AllSamples.trim.contigs.good.unique.good.pick.pick.opti_mcc.0.005.rep.names AllSamples.trim.contigs.good.unique.good.pick.pick.opti_mcc.0.005.rep_modif.names


# Filtering for index-jump and replicates
Rscript Data.treatment.R 

