library(data.table)

# Importing file
entree=read.csv("Tab_final_distri.txt", header=F, sep="\t", dec=".")
colnames(entree)=c("SeqID", "Sample", "Compte")
tableau=entree
tableau=tableau[tableau$Compte!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Replicat=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2]
tableau=as.data.table(tableau)
tableau=tableau[tableau$Sample!="Undetermined",]

# Identification of the maximum number of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg$Compte)
tab_tagJump=tableau[tableau$Compte>=2*maxIndex,]
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, SeqID)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>4,]

# Exporting final contingency table
tab_sortie=dcast(tab_final, SeqID~Sample, value.var="Compte", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_OBITOOLS.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)
