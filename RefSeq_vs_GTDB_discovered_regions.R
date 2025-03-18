setwd("/home/nikolay/WABI/T_van_der_Valk/misc/results_no_sliding_window_no_multimappers")


RefSeq<-read.table(gzfile("757_mammals_RefSeq_BED_microb_contam.txt.gz"),header=TRUE,sep="\t")
GTDB<-read.table(gzfile("757_mammals_GTDB_BED_microb_contam.txt.gz"),header=TRUE,sep="\t")


RefSeq_df<-data.frame(table(RefSeq$ORGANISM_ID))
GTDB_df<-data.frame(table(GTDB$ORGANISM_ID))

intersect_names<-intersect(as.character(RefSeq_df$Var1),as.character(GTDB_df$Var1))
standard_genomes<-read.table("/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/RESULTS_MAMMALS/mammals_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
intersect_names<-intersect_names[intersect_names%in%as.character(standard_genomes$IDS)]

RefSeq_df<-RefSeq_df[as.character(RefSeq_df$Var1)%in%intersect_names,]
GTDB_df<-GTDB_df[as.character(GTDB_df$Var1)%in%intersect_names,]


plot(log10(GTDB_df$Freq+1)~log10(RefSeq_df$Freq+1),
     xlab="log10 ( Number of microbial-like regions discovered with RefSeq pseudo-reads + 1 )",
     ylab="log10 ( Number of microbial-like regions discovered with GTDB pseudo-reads + 1 )",
     #main="Number of discovered regions for mammalian reference genomes with RefSeq and GTDB microbial pseudo-reads",
     main="",
     col="blue",pch=19,cex.axis=1.5,cex.lab=1.5)
abline(coef=c(0,1),col="red",lwd=3)

