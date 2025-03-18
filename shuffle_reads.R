for(k in 1:100000)
{
print(paste0("######################################################## RUNNING ITERATION ",k," ########################################################"))

reads_coords<-read.delim("cr9_67.coords_aligned_reads.bed",header=FALSE,sep="\t")
reads_lengths<-reads_coords$V3-reads_coords$V2
contigs<-read.delim("genome_39321.fna.fai",header=FALSE,sep="\t")

shuffled_read_ids<-vector(); shuffled_read_start<-vector(); shuffled_read_end<-vector()
for(i in 1:length(reads_lengths))
{
random_contig<-sample(as.character(contigs$V1),1)
length_of_random_contig<-contigs[as.character(contigs$V1)==random_contig,]$V2

random_read_start_coord<-sample(1:(length_of_random_contig-reads_lengths[i]), 1)
random_read_end_coord<-random_read_start_coord+reads_lengths[i]

shuffled_read_ids<-append(shuffled_read_ids,random_contig)
shuffled_read_start<-append(shuffled_read_start,random_read_start_coord)
shuffled_read_end<-append(shuffled_read_end,random_read_end_coord)

if(i%%1000==0){print(paste0("Generated ",i," reads"))}
}

shuffled_reads<-data.frame(V1=shuffled_read_ids,V2=shuffled_read_start,V3=shuffled_read_end)

print(head(shuffled_reads))
print(head(reads_coords))

write.table(shuffled_reads,file=paste0("cr9_67.coords_shuffled_reads",k,".bed"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

system("ml bioinfo-tools samtools BEDTools")
system(paste0("sort -k1,1 -k2,2n cr9_67.coords_shuffled_reads",k,".bed > cr9_67.coords_shuffled_reads",k,".sorted.bed"))
system(paste0("bedtools closest -a cr9_67.coords_shuffled_reads",k,".sorted.bed -b coords_micr_contam_GTDB_genome_39321.fna.gz_with_multimappers.sorted.bed -d > cr9_67.coords_shuffled_reads",k,"_annotated_with_closest_contam_region.bed"))
system(paste0("cut -f7 cr9_67.coords_shuffled_reads",k,"_annotated_with_closest_contam_region.bed | awk '{if($1==0)print $0}' | wc -l >> number_of_random_intersects.txt"))
system(paste0("rm cr9_67.coords_shuffled_reads",k,".bed"))
}

