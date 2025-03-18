setwd("/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/")

library("stringr")
library("RColorBrewer")
library("ggplot2")
library("rphylopic")

#pdf("cont_barplot.pdf",paper = "a4",width=210, height=297)

par(mfrow=c(3,2))
par(mar=c(3.8,21,2,2.1))
#par(mar=c(5.1,18,4.1,2.1), mai = c(0.1,2,0.1,1))

mammals<-read.delim("RESULTS_MAMMALS/mammals_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
barplot(rev(mammals$BOC[1:10]*100),names=gsub("_old","",rev(mammals$ORGANISM[1:10])),las=2,
        cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8,horiz=TRUE,col=brewer.pal(6,"Dark2")[1],
        main="MAMMALS",font=3)

uuid <- get_uuid(name = "Rangifer tarandus")
img <- get_phylopic(uuid = uuid)
add_phylopic_base(img = img, x = 1, y = 4.25, height = 7.25)


plants<-read.delim("RESULTS_PLANTS/plants_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
plants_names<-gsub("_"," ",word(word(plants$IDS, 1,2, sep="_"),1,1,sep="-"))
plants_ids<-word(word(plants$IDS,3,3,sep="-"),1,2,sep="_")
plants$REF_IDS<-plants_ids
plants$ORGANISM<-plants_names
#write.table(plants,file="RESULTS_PLANTS/SupplTable2.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
barplot(rev(plants$BOC[1:10]*100),names=rev(plants_names[1:10]),las=2,
        cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8,horiz=TRUE,
        col=brewer.pal(6,"Dark2")[2],main="NCBI REFSEQ PLANTS",font=3)

uuid <- get_uuid(name = "Betula pendula")
img <- get_phylopic(uuid = uuid)
add_phylopic_base(img = img, x = 1.2, y = 3.25, height = 4.25)


vertebrates<-read.delim("RESULTS_VERTEBRATES/verterbrates_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
vertebrates_names<-gsub("_"," ",word(word(vertebrates$IDS, 1,2, sep="_"),1,1,sep="-"))
vertebrates_ids<-word(word(vertebrates$IDS,3,3,sep="-"),1,2,sep="_")
vertebrates$REF_IDS<-vertebrates_ids
vertebrates$ORGANISM<-vertebrates_names
#write.table(vertebrates,file="RESULTS_VERTEBRATES/SupplTable3.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
barplot(rev(vertebrates$BOC[1:10]*100),names=rev(vertebrates_names[1:10]),las=2,
        cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8,horiz=TRUE,col=brewer.pal(6,"Dark2")[3],
        main="VERTEBRATES",xlim=c(0,max(rev(vertebrates$BOC[1:10]*100))+0.04),font=3)

uuid <- get_uuid(name = "Salmo trutta")
img <- get_phylopic(uuid = uuid)
add_phylopic_base(img = img, x = 0.15, y = 3.25, height = 2)


invertebrates<-read.delim("RESULTS_INVERTEBRATES/inverterbrates_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
invertebrates_names<-gsub("_"," ",word(word(invertebrates$IDS, 1,2, sep="_"),1,1,sep="-"))
invertebrates_ids<-word(word(invertebrates$IDS,3,3,sep="-"),1,2,sep="_")
invertebrates$REF_IDS<-invertebrates_ids
invertebrates$ORGANISM<-invertebrates_names
#write.table(invertebrates,file="RESULTS_INVERTEBRATES/SupplTable4.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
barplot(rev(invertebrates$BOC[1:10]*100),names=rev(invertebrates_names[1:10]),las=2,
        cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8,horiz=TRUE,col=brewer.pal(6,"Dark2")[4],
        main="INVERTEBRATES",xlim=c(0,max(rev(invertebrates$BOC[1:10]*100))+0.5),font=3)

uuid <- get_uuid(name = "Drosophila melanogaster")
img <- get_phylopic(uuid = uuid)
add_phylopic_base(img = img, x = 1.5, y = 4.25, height = 3)


phylonorway<-read.delim("RESULTS_PHYLONORWAY/phylonorway_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
phylonorway_annot<-read.delim("RESULTS_PHYLONORWAY/1323_plant_phylonorway_sorted_by_BOC.txt",header=TRUE,sep="\t")
phylonorway$ORGANISM<-phylonorway_annot$ORGANISM[match(as.character(phylonorway$IDS),as.character(phylonorway_annot$IDS))]
phylonorway_names<-phylonorway$ORGANISM
phylonorway$IDS<-gsub(".fna.gz","",gsub("genome_","",as.character(phylonorway$IDS)))
#write.table(phylonorway,file="RESULTS_PHYLONORWAY/SupplTable5.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
barplot(rev(phylonorway$BOC[1:10]*100),names=rev(phylonorway_names[1:10]),las=2,
        xlab="PERCENT OF MICROBIAL CONTAMINATION",
        cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8,
        horiz=TRUE,col=brewer.pal(6,"Dark2")[6],
        main="PHYLONORWAY PLANTS",xlim=c(0,max(rev(phylonorway$BOC[1:10]*100))+1),font=3)

uuid <- get_uuid(name = "Picea abies")
img <- get_phylopic(uuid = uuid)
add_phylopic_base(img = img, x = 50, y = 4.25, height = 8.25)



arthropods<-read.delim("RESULTS_ARTHROPODS/arthropods_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
arthropods_ids<-word(arthropods$IDS,1,2,sep="_")
arthropods$REF_IDS<-arthropods_ids
#write.table(arthropods,file="RESULTS_ARTHROPODS/SupplTable6.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
barplot(rev(arthropods$BOC[1:10]*100),names=rev(arthropods$ORGANISM[1:10]),las=2,
        xlab="PERCENT OF MICROBIAL CONTAMINATION",
        cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8,horiz=TRUE,col=brewer.pal(5,"Dark2")[5],
        main="ARTHROPODS",xlim=c(0,max(rev(arthropods$BOC[1:10]*100))+1),font=3)

uuid <- get_uuid(name = "Latrodectus geometricus")
img <- get_phylopic(uuid = uuid)
add_phylopic_base(img = img, x = 5.5, y = 5, height = 8.25)

#dev.off()

# computing overlapping invertebrates and arthropods
length(intersect(invertebrates_ids,arthropods_ids))#55
overlapping_ids<-intersect(invertebrates_ids,arthropods_ids)
overlapping_invertebrates<-invertebrates[match(intersect(invertebrates_ids,arthropods_ids),as.character(invertebrates$REF_IDS)),]
overlapping_arthropods<-arthropods[match(intersect(invertebrates_ids,arthropods_ids),as.character(arthropods$REF_IDS)),]
head(overlapping_invertebrates)
head(overlapping_arthropods)
plot(overlapping_invertebrates$BOC~overlapping_arthropods$BOC)
abline(c(0,0),c(1,1))

arhropods_new<-arthropods[!as.character(arthropods$REF_IDS)%in%as.character(overlapping_ids),]
write.table(arhropods_new,file="RESULTS_ARTHROPODS/SupplTable6.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
