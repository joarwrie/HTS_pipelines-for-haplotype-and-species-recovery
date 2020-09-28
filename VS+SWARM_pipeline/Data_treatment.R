library(data.table)

# Importing file
entree=read.csv("AllSamples_finalTAB.csv", header=T, sep="\t", dec=".")
tableau=melt(entree, variable.name="UV_Sample", value.name="Nb_reads", id.vars="X.OTU.ID")
tableau=tableau[tableau$Nb_reads!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$UV_Sample), "_"))[,1]
tableau$Replicat=paste(do.call(rbind, strsplit(as.character(tableau$UV_Sample), "_"))[,2], do.call(rbind, strsplit(as.character(tableau$UV_Sample), "_"))[,3], sep="_")
OTU_table=read.csv("AllSamples_OTUlist_modif.txt", header=F, sep="\t", dec=".", stringsAsFactors=F)
colnames(OTU_table)=c("Seed", "SeqID")
tableau=merge(tableau, OTU_table, by.x="X.OTU.ID", by.y="SeqID", all=T)
tableau=as.data.table(tableau)
tableau=tableau[tableau$Echantillon!="Undetermined",]
tableau[,"Somme":=sum(Nb_reads), by=.(Seed, Echantillon)]
tableau=unique(tableau[,4:7])

# Identification of the maximum number of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg$Somme)

# Filtering for index jump
tab_tagJump=tableau[tableau$Somme>=2*maxIndex,]

# Filtering for replicates
tab_tagJump=tab_tagJump[tab_tagJump$Replicat!="ext1_PCR1" & tab_tagJump$Replicat!="ext1_PCR5" & tab_tagJump$Replicat!="ext2_PCR2" & tab_tagJump$Replicat!="ext2_PCR4" & tab_tagJump$Replicat!="ext3_PCR1" & tab_tagJump$Replicat!="ext3_PCR3" & tab_tagJump$Replicat!="ext4_PCR4" & tab_tagJump$Replicat!="ext4_PCR5",]
tab_tagJump$Type="Ethanol"
tab_tagJump[grep("bulk",tab_tagJump$Echantillon),]$Type="Bulk"
tab_tagJump$Bon=TRUE
tab_tagJump[tab_tagJump$Type=="Bulk" & (tab_tagJump$Replicat=="ext4_PCR1" | tab_tagJump$Replicat=="ext4_PCR2" |tab_tagJump$Replicat=="ext4_PCR3"),]$Bon=FALSE
tab_tagJump[tab_tagJump$Type=="Ethanol" & (tab_tagJump$Replicat=="ext1_PCR2" | tab_tagJump$Replicat=="ext1_PCR3" |tab_tagJump$Replicat=="ext1_PCR4"),]$Bon=FALSE
tab_tagJump=tab_tagJump[tab_tagJump$Bon==TRUE,]
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, Seed)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>4,]

# Exporting final contingency table
tab_sortie=dcast(tab_final, Seed~Echantillon, value.var="Somme", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_OBISwarm.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)
