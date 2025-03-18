library("ggridges")
library("ggplot2")
library("ggbreak")
library("scales")
library("ggtext")

setwd("/home/nikolay/WABI/T_van_der_Valk/contaminated_genomes/RESULTS_WITH_MULTIMAPPERS/")

mammals<-read.delim("RESULTS_MAMMALS/mammals_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
plants<-read.delim("RESULTS_PLANTS/plants_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
vertebrates<-read.delim("RESULTS_VERTEBRATES/verterbrates_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
invertebrates<-read.delim("RESULTS_INVERTEBRATES/inverterbrates_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
phylonorway<-read.delim("RESULTS_PHYLONORWAY/phylonorway_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")
arthropods<-read.delim("RESULTS_ARTHROPODS/arthropods_micr_cont_boc_sorted.txt",header=TRUE,sep="\t")

mammals_df<-data.frame(ORGANISM="MAMMALS",BOC=mammals$BOC)
plants_df<-data.frame(ORGANISM="PLANTS",BOC=plants$BOC)
vertebrates_df<-data.frame(ORGANISM="VERTEBRATES",BOC=vertebrates$BOC)
invertebrates_df<-data.frame(ORGANISM="INVERTEBRATES",BOC=invertebrates$BOC)
phylonorway_df<-data.frame(ORGANISM="PHYLONORWAY",BOC=phylonorway$BOC)
arthropods_df<-data.frame(ORGANISM="ARTHROPODS",BOC=arthropods$BOC)

myboc<-rbind(rbind(rbind(rbind(rbind(mammals_df,plants_df),
                         vertebrates_df),invertebrates_df),phylonorway_df),arthropods_df)

myboc$ORGANISM<-factor(myboc$ORGANISM, levels = c("PHYLONORWAY", "ARTHROPODS", "INVERTEBRATES", 
                       "VERTEBRATES", "PLANTS", "MAMMALS"))

p<-ggplot(myboc) +
  geom_density_ridges(data = myboc, aes(x = BOC, y = ORGANISM, fill = ORGANISM)) +
  theme_ridges() + 
  theme(legend.position = "none", 
        axis.title.x=element_textbox_simple(halign=0.5, size=24)) + 
  scale_x_continuous(trans="log10",breaks=trans_breaks("log10", function(x) 10^x)) + 
  xlab("Microbial-like breadth of coverage") + ylab("") + 
  scale_fill_brewer(palette = "Dark2", direction = -1)
#+ geom_phylopic(aes(image = img))
#+ geom_image(aes(image = img))
#add_phylopic(img = img, x = 0.1, y = 0.1)
#+ scale_x_break(c(0.05, 0.6))
p

library("rphylopic")
#library("ggimage")
uuid <- get_uuid(name = "Rangifer tarandus")
img <- get_phylopic(uuid = uuid)
#p + add_phylopic(uuid = uuid, x = 1, y = 1, height = 0.05)

ggplot() +
  scale_x_continuous(trans="log10",breaks=trans_breaks("log10", function(x) 10^x)) +
  add_phylopic(img = img, x = 1.25, y = 1.25, height = 0.05)

