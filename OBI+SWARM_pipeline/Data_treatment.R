library(data.table)

# Importing file
entree=read.csv("Tab_final_distri.txt", header=F, sep="\t", dec=".")
colnames(entree)=c("SeqID", "Sample", "Compte")
OTU_table=read.csv("AllSamples_OTUlist.txt", header=F, sep=" ", dec=".", stringsAsFactors=F)
OTU_table$OTU.ID=paste("OTU", seq(1,nrow(OTU_table),1), sep="_")
OTU_table=melt(OTU_table, variable.name="Autre", value.name="SeqID", id.vars="OTU.ID")
OTU_table$SeqID=gsub(";size=[0-9]+;", "", OTU_table$SeqID)
tableau=merge(entree, OTU_table, by="SeqID", all=T)
tableau=tableau[tableau$Compte!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Replicat=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2]
tableau=as.data.table(tableau)
tableau=tableau[tableau$Sample!="Undetermined",]
tableau[,"Somme":=sum(Compte), by=.(OTU.ID, Sample)]
tableau=unique(tableau[,c(2,4,6,7,8)])

# Identification of the maximum number of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg$Somme)
tab_tagJump=tableau[tableau$Somme>=2*maxIndex,]
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, OTU.ID)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>4,]

# Exporting final contingency table
tab_sortie=dcast(tab_final, OTU.ID~Sample, value.var="Somme", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_OBISwarm.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)
