library(data.table)

# Filtering for index-jumping and replicates
entree=read.csv("AllSamples_finalTAB.csv", header=T, sep="\t", dec=".")
tableau=melt(entree, variable.name="UV_Sample", value.name="Nb_reads", id.vars="X.OTU.ID")
tableau=tableau[tableau$Nb_reads!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$UV_Sample), "_"))[,1]
tableau$Replicat=paste(do.call(rbind, strsplit(as.character(tableau$UV_Sample), "_"))[,2], do.call(rbind, strsplit(as.character(tableau$UV_Sample), "_"))[,3], sep="_")
tableau=as.data.table(tableau)
tableau=tableau[tableau$Echantillon!="Undetermined",]
tableau[,"Compte":=sum(Nb_reads), by=.(X.OTU.ID, Echantillon, Replicat)]
tableau=unique(tableau[,c(1,4,5,6)])

# Identification of the maximum number of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg$Compte)
tab_tagJump=tableau[tableau$Compte>=maxIndex,]
tab_tagJump=tab_tagJump[tab_tagJump$Replicat!="ext1_PCR1" & tab_tagJump$Replicat!="ext1_PCR5" & tab_tagJump$Replicat!="ext2_PCR2" & tab_tagJump$Replicat!="ext2_PCR4" & tab_tagJump$Replicat!="ext3_PCR1" & tab_tagJump$Replicat!="ext3_PCR3" & tab_tagJump$Replicat!="ext4_PCR4" & tab_tagJump$Replicat!="ext4_PCR5",]
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, X.OTU.ID)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>4,]
tab_final$Sample=paste(tab_final$Echantillon, tab_final$Replicat, sep="_")

# Exporting final contingency table
tab_sortie=dcast(tab_final, X.OTU.ID~Sample, value.var="Compte", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_VSSwarm.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)
