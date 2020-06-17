##################
# PIPELINE DADA2
##################

# Require R-3.5.1
# Require dada2-1.13.1

library(dada2)
library(data.table)

# fastq files reading
path <- "./Demult_files"

fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names=T))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

# Quality visualisation
plotQualityProfile(fnFs[1:2])

# Filtering and cutting reads according to their quality
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(262,170), maxN=0, truncQ=0, multithread=T, rm.phix=F)

# Error modelling
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)

# Data dereplication
derepF <- derepFastq(filtFs)
derepR <- derepFastq(filtRs)

# Denoising
dadaF <- dada(derepF, err=errF, multithread=T, pool=T)
dadaR <- dada(derepR, err=errR, multithread=T, pool=T)

# Merging read pairs
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap=20, maxMismatch=0)

# Length selection
seqtab <- makeSequenceTable(mergers)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 402:422]

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=T, verbose=T)

# Exporting data files
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(file="tab_track_pool_dada2.csv", track, col.names=T, row.names=T, sep="\t", dec=".", quote=F)

ASVs=data.frame(Code=paste(rep("ASV",length(colnames(seqtab2))), 1:length(colnames(seqtab2)), sep="_"), Sequence=getSequences(seqtab2))
ASVs.nochim=data.frame(Code=paste(rep("ASV",length(colnames(seqtab.nochim))), 1:length(colnames(seqtab.nochim)), sep="_"), Sequence=getSequences(seqtab.nochim))
colnames(seqtab2)=paste(rep("ASV",length(colnames(seqtab2))), 1:length(colnames(seqtab2)), sep="_")
colnames(seqtab.nochim)=paste(rep("ASV",length(colnames(seqtab.nochim))), 1:length(colnames(seqtab.nochim)), sep="_")
write.table(file="List_ASVs_withChim.fasta", ASVs, col.names=F, row.names=F, sep="\n", dec=".", quote=F)
write.table(file="tab_distri_ASVs_withChim.csv", seqtab2, col.names=T, row.names=T, sep="\t", dec=".", quote=F)
write.table(file="List_ASVs_noChim.fasta", ASVs.nochim, col.names=F, row.names=F, sep="\n", dec=".", quote=F)
write.table(file="tab_distri_ASVs_noChim.csv", seqtab.nochim, col.names=T, row.names=T, sep="\t", dec=".", quote=F)

# Contingency table modification
entree=as.data.frame(seqtab.nochim)
entree$Sample=row.names(entree)
tableau=melt(entree, variable.name="ASV", value.name="Compte", id.vars="Sample")
tableau=tableau[tableau$Compte!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Replicat=paste(do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2], do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,3], sep="_")
tableau=as.data.table(tableau)
tableau=tableau[tableau$Sample!="Undetermined",]

# Identification of the maximum number of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg$Compte)
tab_tagJump=tableau[tableau$Compte>=maxIndex,]
tab_tagJump=tab_tagJump[tab_tagJump$Replicat!="ext1_PCR1" & tab_tagJump$Replicat!="ext1_PCR5" & tab_tagJump$Replicat!="ext2_PCR2" & tab_tagJump$Replicat!="ext2_PCR4" & tab_tagJump$Replicat!="ext3_PCR1" & tab_tagJump$Replicat!="ext3_PCR3" & tab_tagJump$Replicat!="ext4_PCR4" & tab_tagJump$Replicat!="ext4_PCR5",]
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, ASV)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>4,]

# Exporting final contingency table
tab_sortie=dcast(tab_final, ASV~Sample, value.var="Compte", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_DADA2.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)
