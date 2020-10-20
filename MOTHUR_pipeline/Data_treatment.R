library(data.table)

# Importing file
entree=read.csv("AllSamples_final_table.shared", header=T, sep="\t", dec=".", stringsAsFactors = F)
tableau=melt(entree, variable.name="Sample", value.name="Nb_reads", id.vars=c("Representative_Sequence", "OTU", "total"))
tableau=tableau[tableau$Nb_reads!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Replicat=gsub(".COIBotr", "", do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2])

OTU_table=read.csv("AllSamples.trim.contigs.good.unique.good.pick.pick.opti_mcc.0.005.rep_modif.names", header=F, sep="\t", dec=".")
colnames(OTU_table)=c("Seed", "SeqID")
tableau=merge(tableau, OTU_table, by.x="Nouvelle", by.y="SeqID", all=T)
tableau=as.data.table(tableau)
tableau=tableau[tableau$Echantillon!="Undetermined",]
tableau[,"Somme":=sum(Nb_reads), by=.(Seed, Sample)]
tableau=unique(tableau[,6:10])

# Identification of the maximum number of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg$Somme)

# Filtering for index jump
tab_tagJump=tableau[tableau$Somme>=2*maxIndex,]

# Filtering for replicates
tab_tagJump=unique(tab_tagJump[,2:5])
tab_tagJump$Type="Ethanol"
tab_tagJump[grep("bulk",tab_tagJump$Echantillon),]$Type="Bulk"
tab_tagJump$Bon=TRUE
tab_tagJump[tab_tagJump$Type=="Bulk" & (tab_tagJump$Replicat=="rep41" | tab_tagJump$Replicat=="rep42" |tab_tagJump$Replicat=="rep43"),]$Bon=FALSE
tab_tagJump[tab_tagJump$Type=="Ethanol" & (tab_tagJump$Replicat=="rep12" | tab_tagJump$Replicat=="rep13" |tab_tagJump$Replicat=="rep14"),]$Bon=FALSE
tab_tagJump=tab_tagJump[tab_tagJump$Bon==TRUE,]
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, Seed)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>4,]

# Exporting final contingency table
tab_sortie=dcast(tab_final, Seed~Echantillon, value.var="Somme", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_mothur.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)