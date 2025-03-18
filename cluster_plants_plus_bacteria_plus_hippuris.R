setwd("/home/nikolay/WABI/T_van_der_Valk/IGV_PhyloNorway/")

#Read distance matrix and shorten reference genome ids
df<-read.delim("dist_plants_plus_bacteria_plus_hippuris.txt",
               header=TRUE,row.names=1,sep="\t")
df<-df[-which(grepl("Bathycoccus_prasinos",rownames(df))),]; 
df[,"Bathycoccus_prasinos_strain_RCC_1105.tax41875.GCF_002220235.1_ASM222023v1_genomic.dustmasked.fna"]<-NULL
df<-df[-which(grepl("Ostreococcus_tauri",rownames(df))),]; 
df[,"Ostreococcus_tauri_strain_RCC4221.tax70448.GCF_000214015.3_version_140606_genomic.dustmasked.fna"]<-NULL
df<-df[-which(grepl("Chlamydomonas_reinhardtii",rownames(df))),]; 
df[,"Chlamydomonas_reinhardtii_strain_CC_503_cw92_mt.tax3055.GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.dustmasked.fna"]<-NULL

library("stringr")
plants_names<-word(word(colnames(df),1,1,sep=".tax"), 1,2, sep="_")
plants_names[grepl("^GCF_",plants_names)]<-tolower(plants_names[grepl("^GCF_",plants_names)])
plants_names[grepl("contam",plants_names)]<-"Hippuris_vulgaris_microbial_like"
plants_names[grepl("clean",plants_names)]<-"Hippuris_vulgaris_endogenous"
colnames(df)<-plants_names; rownames(df)<-plants_names
df[1:5,1:5]

color_vec<-ifelse(grepl("^gcf_",colnames(df)),"darkred","darkgreen")
color_vec[grepl("Hippuris_vulgaris_endogenous",colnames(df))]<-"magenta"
color_vec[grepl("Hippuris_vulgaris_microbial_like",colnames(df))]<-"darkorange"

#Run hierarchical clustering and plot dendrogram
library("ggplot2"); library("ggdendro")
#model<-hclust(as.dist(df),"ward.D2")
PC <- prcomp(t(log10(df + 1)), center=FALSE, scale=FALSE)
model<-hclust(dist(PC$x[,1:2]),"ward.D2")
dhc <- as.dendrogram(model)
nodePar <- list(lab.cex = 0.5, pch = c(NA, 19), cex = 0.5, col = "blue")
library("dendextend")
color_vec <- color_vec[order.dendrogram(dhc)]
labels_colors(dhc)<-color_vec
pdf("dendrogram_portrate.pdf",paper = "a4",width=210, height=297)
par(oma=c(0,0,0,0), mar=c(0,0,0,10))
par(cex=0.5)
plot(dhc, ylab = "", nodePar = nodePar, horiz = TRUE)
legend("topleft", inset=.02,
       #c("Bacterial RefSeq genomes","Plant RefSeq genomes","Hippuris vulgaris endogenous regions","Hippuris vulgaris microbial-like regions"),
       c("Bacterial RefSeq genomes","Plant RefSeq genomes",
         expression(paste(italic("Hippuris vulgaris")," endogenous regions")),
         expression(paste(italic("Hippuris vulgaris")," microbial-like regions"))),
       fill=c("darkred","darkgreen","magenta","darkorange"), horiz=FALSE, cex=2)
dev.off()



# Run PCA
PC <- prcomp(t(log10(df + 1)), center=FALSE, scale=FALSE)
expl_var <- PC$sdev^2/sum(PC$sdev^2)
barplot(expl_var[1:20],ylab="EXPLAINED VARIANCE",
        main="VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS",
        names.arg=paste0("PC",seq(1:20)),col="darkgreen")
plot(PC$x[,1:2], main=paste0("PCA PLOT: PLANTS + BACTERIA ORGANISMS"), 
     xlab="PC1", ylab="PC2", type="n")
text(PC$x[,1:2], labels=paste(colnames(df),sep=""),col=color_vec,cex=0.5)

my_mds <- cmdscale(df)
plot(my_mds, main=paste0("MDS PLOT: PLANTS + BACTERIA ORGANISMS"), 
     xlab="MDS1", ylab="MDS2", type="n")
text(my_mds, labels=paste(colnames(df),sep=""),col=color_vec,cex=0.5)

# select informative PCs
print("SELECTING SIGNIFICANT PRINCIPAL COMPONENTS")
N_perm<-10
expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = N_perm)
for(k in 1:N_perm)
{
  df_perm <- apply(df,2,sample)
  PC_perm <- prcomp(t(log10(df_perm+1)), center=TRUE, scale=FALSE)
  expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
  print(paste0("FINISHED ",k," PERMUTATIONS"))
}
plot(expl_var[1:20]~seq(1:20),ylab="EXPLAINED VARIANCE",
     main="VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS",
     col="darkgreen",type='o',xlab="PRINCIPAL COMPONENTS")
lines(colMeans(expl_var_perm)[1:20]~seq(1:20),col="red")
legend("topright",c("Explained by PCS","Explained by chance"),fill=c("darkgreen","red"),
       inset=0.02)
pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
plot(pval[1:20]~seq(1:20),col="darkred",type='o',xlab="PRINCIPAL COMPONENTS",ylab="PVALUE",
     main="SIGNIFICANCE OF PRINCIPAL COMPONENTS")
optPC<-head(which(pval>=0.05),1)-1
mtext(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))
print(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC,
             ", THEY TOGETHER EXPLAIN ",round(sum(expl_var[1:optPC])*100,0),"% OF VARIANCE"))

############################################# tSNE ANALYSIS ################################
library("Rtsne")
optPerp<-round(sqrt(dim(df)[2]),0)
tsne_opt_perp <- Rtsne(t(log10(df + 1)),initial_dims=optPC,verbose=TRUE,
                       check_duplicates=FALSE,perplexity=optPerp,dims=2,max_iter=10000)
plot(tsne_opt_perp$Y,main=paste0("TSNE PLOT: PLANTS + BACTERIA"),
     xlab="tSNE1",ylab="tSNE2", type="n")
text(tsne_opt_perp$Y,labels=paste(colnames(df),sep=""),col=color_vec,cex=0.5)

