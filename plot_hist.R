setwd("/home/nikolay/WABI/T_van_der_Valk/Hippuris/with_multimappers/")

aligned_read<-read.delim("cr9_67.coords_aligned_reads.sorted.bed",header=FALSE,sep="\t")
head(aligned_read)

aligned_reads_KapK<-
  read.delim("KapK/with_multimappers/69_B2_100_L0_KapK-12-1-35_Ext-12_Lib-12.coords_aligned_reads.sorted.bed",
                               header=FALSE,sep="\t")
head(aligned_reads_KapK)

fai<-read.delim("cont_fraction_per_Hippuris_contig.txt",header=TRUE,sep="\t")
head(fai)

#par(mfrow = c(1,2))
hist(fai$CONT_FRACTION,breaks=50,col="darkgreen",
     #xlab="Contamination fraction of Hippuris vulgaris contig",
     xlab=expression(paste("Microbial-like fraction of ",italic("Hippuris vulgaris")," contig")),
     #main="All 433 631 Hippuris vulgaris contigs",
     main="",
     cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8)

fai_with_aligned_reads<-fai[as.character(fai$CONTIG)%in%unique(as.character(aligned_read$V1)),]

hist(fai_with_aligned_reads$CONT_FRACTION,breaks=40,col="darkred",
     #xlab="Microbial-like fraction of Hippuris vulgaris contig",
     xlab=expression(paste("Microbial-like fraction of ",italic("Hippuris vulgaris")," contig")),
     #main="Only 20 213 Hippuris vulgaris contigs where cr9_67 sequenced reads align",
     main="",
     cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8)

text(x = 0, y = 17500, labels = "A)", xpd = NA, cex=3)


fai_with_aligned_reads_KapK<-fai[as.character(fai$CONTIG)%in%unique(as.character(aligned_reads_KapK$V1)),]
hist(fai_with_aligned_reads_KapK$CONT_FRACTION,breaks=40,col="darkred",
     #xlab="Microbial-like fraction of Hippuris vulgaris contig",
     xlab=expression(paste("Microbial-like fraction of ",italic("Hippuris vulgaris")," contig")),
     #main="Only 73 911 Hippuris vulgaris contigs where 69_B2_100_L0_KapK sequenced reads align",
     main="",
     cex.names=1.8,cex.axis=1.5,cex.lab=1.5,cex.main=1.8)

text(x = 0, y = 58000, labels = "B)", xpd = NA, cex=3)



# Random read assignemnts and intersects with microbial-like regions in Hippuris vulgaris
intersect_numbers<-as.numeric(scan("number_of_random_intersects.txt",what="character"))
hist(intersect_numbers/119854,xlab="Fraction of reads",
     #main="Arctic sample cr9_67 from Wang et al., Nature 2021",
     main="",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.8, breaks = 30, 
     col = "darkorange")
#mtext("Histogram of randomly placed reads intersecting contaminated regions in Hippuris vulgaris reference")
text(x = 0.585, y = 25, labels = "A)", xpd = NA, cex=3)


intersect_numbers<-as.numeric(scan("KapK/with_multimappers/number_of_random_intersects.txt",
                                   what="character"))
hist(intersect_numbers/1367627,xlab="Fraction of reads",
     main="Greenland sample 69_B2_100_L0_KapK-12-1-35 from Kjajer et al., Nature 2022",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.8, col = "darkmagenta")
mtext("Histogram of randomly placed reads intersecting contaminated regions in Hippuris vulgaris reference")


hist(intersect_numbers/1367627,xlab="Fraction of reads",
     #main="Greenland sample 69_B2_100_L0_KapK-12-1-35 from Kjajer et al., Nature 2022",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.8, col = "darkmagenta",add=TRUE)
mtext("Histogram of randomly placed reads intersecting contaminated regions in Hippuris vulgaris reference")


