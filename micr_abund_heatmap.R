##################################### MAMMALS #####################################################
setwd("/home/nikolay/WABI/T_van_der_Valk/Manuscript/contamination_paper/micr_abundance/mammals/")

abund<-read.delim("micr_abund_matrix_top10.txt",header=TRUE,
                  row.names=1,check.names=FALSE,sep="\t")
#rownames(abund)<-abund$MICROBE; abund$MICROBE<-NULL
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

Path<-"/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/"
annot<-read.delim(paste0(Path,"RESULTS_MAMMALS/mammals_micr_cont_boc_sorted.txt"),
                  header=TRUE,sep="\t")
annot<-annot[1:200,]
abund<-subset(abund,select=colnames(abund)[colnames(abund)%in%annot$IDS])
colnames(abund)<-annot$ORGANISM[match(as.character(colnames(abund)),as.character(annot$IDS))]
micr_annot<-read.delim("/home/nikolay/WABI/T_van_der_Valk/scripts/GTDB_fna2name.txt",
                       header=FALSE,sep="\t")
rownames(abund)<-as.character(micr_annot$V3)[match(rownames(abund),as.character(micr_annot$V1))]
abund[1:5,1:5]

#remove rare microbes
abund<-abund[rowSums(abund!=0)>3,]
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])

library("pheatmap")
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)

#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(abund)[2])
{
  abund[,i]<-abund[,i]/sum(abund[,i],na.rm=TRUE)
}
abund[1:5,1:5]
pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(abund,cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(abund,cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)
pdf("micr_abund_heatmap_portrate.pdf",paper = "a4",width=210, height=297)
pheatmap(t(abund),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=5)
dev.off()


#extract colnames and rownames from clustered heatmap
out<-pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8,
              clustering_method="complete",show_rownames=T)
rownames(abund[out$tree_row[["order"]],])
colnames(abund[out$tree_col[["order"]],])

plot(out$tree_col)

###################################################################################################

##################################### PLANTS ######################################################
setwd("/home/nikolay/WABI/T_van_der_Valk/Manuscript/contamination_paper/micr_abundance/plants/")

abund<-read.delim("micr_abund_matrix_top10_plants.txt",header=TRUE,
                  row.names=1,check.names=FALSE,sep="\t")
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

Path<-"/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/"
annot<-read.delim(paste0(Path,"RESULTS_PLANTS/plants_micr_cont_boc_sorted.txt"),
                  header=TRUE,sep="\t")
library("stringr")
colnames(abund)<-gsub("_"," ",word(word(colnames(abund), 1,2, sep="_"),1,1,sep="-"))
micr_annot<-read.delim("/home/nikolay/WABI/T_van_der_Valk/scripts/GTDB_fna2name.txt",
                       header=FALSE,sep="\t")
rownames(abund)<-as.character(micr_annot$V3)[match(rownames(abund),as.character(micr_annot$V1))]
abund[1:5,1:5]

library("pheatmap")
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)

#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(abund)[2])
{
  abund[,i]<-abund[,i]/sum(abund[,i],na.rm=TRUE)
}
abund[1:5,1:5]
pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(abund,cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(abund,cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)
pdf("micr_abund_heatmap_portrate_plants.pdf",paper = "a4",width=210, height=297)
pheatmap(t(abund),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=5)
dev.off()


#extract colnames and rownames from clustered heatmap
out<-pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8,
              clustering_method="complete",show_rownames=T)
rownames(abund[out$tree_row[["order"]],])
colnames(abund[out$tree_col[["order"]],])

plot(out$tree_col)

###############################################################################################

##################################### INVERTEBRATES ###########################################
setwd("/home/nikolay/WABI/T_van_der_Valk/Manuscript/contamination_paper/micr_abundance/invertebrates/")

abund<-read.delim("micr_abund_matrix_top10_invertebrates.txt",header=TRUE,
                  row.names=1,check.names=FALSE,sep="\t")
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

Path<-"/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/"
annot<-read.delim(paste0(Path,"RESULTS_INVERTEBRATES/inverterbrates_micr_cont_boc_sorted.txt"),
                  header=TRUE,sep="\t")
library("stringr")
colnames(abund)<-gsub("_"," ",word(word(colnames(abund), 1,2, sep="_"),1,1,sep="-"))
micr_annot<-read.delim("/home/nikolay/WABI/T_van_der_Valk/scripts/GTDB_fna2name.txt",
                       header=FALSE,sep="\t")
rownames(abund)<-as.character(micr_annot$V3)[match(rownames(abund),as.character(micr_annot$V1))]
abund[1:5,1:5]

library("pheatmap")
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)

#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(abund)[2])
{
  abund[,i]<-abund[,i]/sum(abund[,i],na.rm=TRUE)
}
abund[1:5,1:5]
pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=7)
pheatmap(abund,cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(abund,cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)
pdf("micr_abund_heatmap_portrate_invertebrates.pdf",paper = "a4",width=210, height=297)
pheatmap(t(abund),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=4)
dev.off()


#extract colnames and rownames from clustered heatmap
out<-pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8,
              clustering_method="complete",show_rownames=T)
rownames(abund[out$tree_row[["order"]],])
colnames(abund[out$tree_col[["order"]],])

plot(out$tree_col)

################################################################################################

####################################### VERTEBRATES ###########################################
setwd("/home/nikolay/WABI/T_van_der_Valk/Manuscript/contamination_paper/micr_abundance/vertebrates/")

abund<-read.delim("micr_abund_matrix_top10_vertebrates.txt",header=TRUE,
                  row.names=1,check.names=FALSE,sep="\t")
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

Path<-"/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/"
annot<-read.delim(paste0(Path,"RESULTS_VERTEBRATES/verterbrates_micr_cont_boc_sorted.txt"),
                  header=TRUE,sep="\t")
library("stringr")
colnames(abund)<-gsub("_"," ",word(word(colnames(abund), 1,2, sep="_"),1,1,sep="-"))
micr_annot<-read.delim("/home/nikolay/WABI/T_van_der_Valk/scripts/GTDB_fna2name.txt",
                       header=FALSE,sep="\t")
rownames(abund)<-as.character(micr_annot$V3)[match(rownames(abund),as.character(micr_annot$V1))]
abund[1:5,1:5]

library("pheatmap")
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)

#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(abund)[2])
{
  abund[,i]<-abund[,i]/sum(abund[,i],na.rm=TRUE)
}
abund[1:5,1:5]
pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=7)
pheatmap(abund,cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(abund,cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)
pdf("micr_abund_heatmap_portrate_vertebrates.pdf",paper = "a4",width=210, height=297)
pheatmap(t(abund),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=4)
dev.off()


#extract colnames and rownames from clustered heatmap
out<-pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8,
              clustering_method="complete",show_rownames=T)
rownames(abund[out$tree_row[["order"]],])
colnames(abund[out$tree_col[["order"]],])

plot(out$tree_col)

#################################################################################################

####################################### ARTHROPODS ###########################################
setwd("/home/nikolay/WABI/T_van_der_Valk/Manuscript/contamination_paper/micr_abundance/arthropods/")

abund<-read.delim("micr_abund_matrix_top10_arthropods.txt",header=TRUE,
                  row.names=1,check.names=FALSE,sep="\t")
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

Path<-"/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/"
annot<-read.delim(paste0(Path,"RESULTS_ARTHROPODS/arthropods_micr_cont_boc_sorted.txt"),
                  header=TRUE,sep="\t")
annot<-annot[1:200,]
abund<-subset(abund,select=colnames(abund)[colnames(abund)%in%annot$IDS])
colnames(abund)<-annot$ORGANISM[match(as.character(colnames(abund)),as.character(annot$IDS))]
micr_annot<-read.delim("/home/nikolay/WABI/T_van_der_Valk/scripts/GTDB_fna2name.txt",
                       header=FALSE,sep="\t")
rownames(abund)<-as.character(micr_annot$V3)[match(rownames(abund),as.character(micr_annot$V1))]

#remove rare microbes
abund<-abund[rowSums(abund!=0)>3,]
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

library("pheatmap")
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)

#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(abund)[2])
{
  abund[,i]<-abund[,i]/sum(abund[,i],na.rm=TRUE)
}
abund[1:5,1:5]
pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=7)
pheatmap(abund,cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(abund,cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)
pdf("micr_abund_heatmap_portrate_arthropods.pdf",paper = "a4",width=210, height=297)
pheatmap(t(abund),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=4)
dev.off()


#extract colnames and rownames from clustered heatmap
out<-pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8,
              clustering_method="complete",show_rownames=T)
rownames(abund[out$tree_row[["order"]],])
colnames(abund[out$tree_col[["order"]],])

plot(out$tree_col)

################################################################################################


####################################### PHYLONORWAY ###########################################
setwd("/home/nikolay/WABI/T_van_der_Valk/Manuscript/contamination_paper/micr_abundance/phylonorway/")

abund<-read.delim("micr_abund_matrix_top10_phylonorway.txt",header=TRUE,
                  row.names=1,check.names=FALSE,sep="\t")
abund<-subset(abund,select=colnames(abund)[colSums(abund)!=0])
abund[1:5,1:5]

Path<-"/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/"
annot<-read.delim(paste0(Path,"RESULTS_PHYLONORWAY/1323_plant_phylonorway_sorted_by_BOC.txt"),
                  header=TRUE,sep="\t")
annot<-annot[1:200,]
abund<-subset(abund,select=colnames(abund)[colnames(abund)%in%as.character(annot$IDS)])
colnames(abund)<-as.character(annot$ORGANISM)[match(colnames(abund),as.character(annot$IDS))]
micr_annot<-read.delim("/home/nikolay/WABI/T_van_der_Valk/scripts/GTDB_fna2name.txt",
                       header=FALSE,sep="\t")
rownames(abund)<-as.character(micr_annot$V3)[match(rownames(abund),as.character(micr_annot$V1))]

#remove rare microbes
abund<-abund[,is.na(colnames(abund))==FALSE]
abund<-abund[rowSums(abund!=0)>0,]
abund<-subset(abund,select=colnames(abund)[colSums(abund,na.rm=TRUE)!=0])
abund[1:5,1:5]

library("pheatmap")
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(log10(abund+1),cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)

#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(abund)[2])
{
  abund[,i]<-abund[,i]/sum(abund[,i],na.rm=TRUE)
}
abund[1:5,1:5]
pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=7)
pheatmap(abund,cluster_rows=TRUE,cluster_cols=FALSE,fontsize=8)
pheatmap(abund,cluster_rows=FALSE,cluster_cols=TRUE,fontsize=8)
pdf("micr_abund_heatmap_portrate_phylonorway.pdf",paper = "a4",width=210, height=297)
pheatmap(t(abund),cluster_rows=TRUE,cluster_cols=TRUE,fontsize=4)
dev.off()


#extract colnames and rownames from clustered heatmap
out<-pheatmap(abund,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=8,
              clustering_method="complete",show_rownames=T)
rownames(abund[out$tree_row[["order"]],])
colnames(abund[out$tree_col[["order"]],])

plot(out$tree_col)


