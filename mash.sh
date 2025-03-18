
# get contaminated sequences
bedtools getfasta -fi genome_39321.fna -bed contam_coords_genome_39321.fna.gz.bed -fo contam_seqs_genome_39321.fna

# get clean sequences, extract coordinates of clean regions with bedtools complement
cut -f1,2 genome_39321.fna.fai > genome.fai
bedtools complement -i contam_coords_genome_39321.fna.gz.bed -g genome.fai > clean_coords_genome_39321.fna.gz.bed
bedtools getfasta -fi genome_39321.fna -bed clean_coords_genome_39321.fna.gz.bed -fo clean_seqs_genome_39321.fna


# run mash to compute k-mer distance
ml bioinfo-tools mash

mash sketch clean_seqs_genome_39321.fna
mash sketch contam_seqs_genome_39321.fna
mash dist clean_seqs_genome_39321.fna.msh contam_seqs_genome_39321.fna.msh



# sampling 10 000 sequences from clean and contaminated collection
ml seqtk

for i in {1..100}
do
echo ${i}
seqtk sample -s${i} clean_seqs_genome_39321.fna 10000 > sample_seqs/clean_${i}.fna
seqtk sample -s${i} contam_seqs_genome_39321.fna 10000 > sample_seqs/contam_${i}.fna
done


# sketching all reference genomes (plants, bacteri and Hippuris clean and contaminated sequences) in the foder
for i in $(ls *.fna)
do
echo ${i}
mash sketch ${i}
done


# computing mash k-mer distance between all pairs of referencec genomes
for i in $(ls *.fna)
do
echo ${i}
for j in $(ls *.fna)
do
mash dist ${i}.msh ${j}.msh >> /proj/snic2021-23-584/nobackup/nikolay/IGV_PhyloNorway/mash_dist_plants_plus_bacteria_plus_hippuris.txt
done
done

