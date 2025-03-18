############################################## OPPOSSUM ##############################################
k<-c(0,5,10,25,50)

nreads_RefSeq<-c(584,11441,12835,12965,13090)
nreads_GTDB<-c(5085,52496,62445,63524,63989)

plot(nreads_GTDB~k,xlab="MAXIMUM NUMBER OF MULTIMAPPERS TO ALLOW",ylab="TOTAL NUMBER OF ALIGNED READS",
     type="o",pch=19,col="red",ylim=c(0,70000),main="Monodelphis domestica: number of aligned reads when including multimappers")
lines(nreads_RefSeq~k,col="blue",type="o",pch=19)

legend("topleft", inset=.02, c("RefSeq: 481 + 2400 unique","GTDB: 2014 + 10940 unique"), fill=c("blue","red"))



nregions_GTDB<-c(243,719,710,715,714)
nregions_RefSeq<-RefSeq<-c(92,221,218,215,216)

plot(nregions_GTDB~k,xlab="MAXIMUM NUMBER OF MULTIMAPPERS TO ALLOW",ylab="NUMBER OF DISCOVERED REGIONS WITH MICROBIAL CONTAMINATiON",
     type="o",pch=19,col="red",ylim=c(0,800),main="Monodelphis domestica: number of contaminated regions when including multimappers")
lines(nregions_RefSeq~k,col="blue",type="o",pch=19)

legend("bottomright", inset=.02, c("RefSeq: k5_138, k10_203, k25_213, k50_213 overlap","GTDB: k5_535, k10_635, k25_693, k50_693 overlap"), 
       fill=c("blue","red"))


########################################### LOXODONTA AFRICANA 3 #####################################################
k<-c(0,5,10,25,50)

nreads_RefSeq<-c(920,4031,4149,4336,4504)

plot(nreads_RefSeq~k,xlab="MAXIMUM NUMBER OF MULTIMAPPERS TO ALLOW",ylab="TOTAL NUMBER OF ALIGNED READS",
     type="o",pch=19,col="blue",ylim=c(0,10000),main="Loxodonta Africana 3: number of aligned reads when including multimappers")

legend("topleft", inset=.02, c("RefSeq"), fill=c("blue"))



nregions_RefSeq<-RefSeq<-c(69,238,239,240,241)

plot(nregions_RefSeq~k,xlab="MAXIMUM NUMBER OF MULTIMAPPERS TO ALLOW",ylab="NUMBER OF DISCOVERED REGIONS WITH MICROBIAL CONTAMINATiON",
     type="o",pch=19,col="blue",ylim=c(0,800),main="Loxodonta Africana 3: number of contaminated regions when including multimappers")

legend("topleft", inset=.02, c("RefSeq"), fill=c("blue"))


############################################## COMBINE AND PLOT VIA GGPLOT #############################################

k<-c(0,5,10,25,50)

nreads_opossum<-c(584,11441,12835,12965,13090)
nreads_loxodonta<-c(920,4031,4149,4336,4504)

#par(mfrow=c(2,1))

plot(nreads_opossum~k,xlab="MAXIMUM NUMBER OF MULTI-MAPPERS TO ALLOW",ylab="NUMBER OF ALIGNED MICROBIAL PSEUDO-READS",
     type="o",pch=19,col="red",ylim=c(0,15000),cex.axis=1.3,cex.lab=1.3,cex.main=1.8,lwd=3)

lines(nreads_loxodonta~k,type="o",pch=19,col="blue",lwd=3)

legend("bottomright", inset=.02, c(expression(paste(italic("Loxodonta africana"),", GCF_000001905.1")),
                                   expression(paste(italic("Monodelphis domestica"),", GCA_027887165.1"))), 
       fill=c("blue","red"), cex = 1.5)
text(x = 0.5, y = 14500, labels = "A)", xpd = NA, cex=3)


k<-c(0,5,10,25,50)

nregions_opossum<-c(92,221,218,215,216)
nregions_loxodonta<-c(69,238,239,240,241)

plot(nregions_opossum~k,xlab="MAXIMUM NUMBER OF MULTI-MAPPERS TO ALLOW",ylab="NUMBER OF DETECTED MICROBIAL-LIKE REGIONS",
     type="o",pch=19,col="red",ylim=c(50,300),cex.axis=1.3,cex.lab=1.3,cex.main=1.8,lwd=3)

lines(nregions_loxodonta~k,type="o",pch=19,col="blue",lwd=3)

legend("bottomright", inset=.02, c(expression(paste(italic("Loxodonta africana"),", GCF_000001905.1")),
                                   expression(paste(italic("Monodelphis domestica"),", GCA_027887165.1"))), 
       fill=c("blue","red"), cex = 1.5)
text(x = 1, y = 290, labels = "B)", xpd = NA, cex=3)
