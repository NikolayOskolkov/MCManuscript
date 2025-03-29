# Load necessary library
library("reshape2")

# Read data
df<-read.table(gzfile("/home/nikolay/WABI/T_van_der_Valk/IGV_PhyloNorway/mash_phylonorway_endo_bact.dist.txt.gz"),
               header=FALSE,sep="\t")
colnames(df)<-c('RefID', 'QueryID', 'MashDistance', 'Pvalue', 'MatchHash')
head(df)

# Get all unique feature names
features <- unique(c(df$RefID, df$QueryID))

# Create an empty matrix
distance_matrix <- matrix(0, nrow = length(features), ncol = length(features),
                          dimnames = list(features, features))

# Fill in the distances
for (i in 1:nrow(df)) {
  f1 <- df$RefID[i]
  f2 <- df$QueryID[i]
  dist <- df$MashDistance[i]
  distance_matrix[f1, f2] <- dist
  distance_matrix[f2, f1] <- dist  # Ensure symmetry
}

# Print matrix
distance_matrix[1:5,1:5]
write.table(distance_matrix, file="/home/nikolay/WABI/T_van_der_Valk/IGV_PhyloNorway/dist_phylonorway_endo_bact.txt", 
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep="\t")
