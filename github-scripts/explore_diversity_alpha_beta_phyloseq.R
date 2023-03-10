############ Initial bar plots and alpha/beta diversity ###########
## Dec 2022
## https://joey711.github.io/phyloseq/plot_bar-examples.html


setwd("~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq")

##Installing Qiime2R
#qiime2R is currently available via github which can easily be installed in R via the following command:
##devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")

library(qiime2R)
library(decontam)
library(biomformat);packageVersion("biomformat")
library(ggplot2)
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library(microViz)
library(vegan)
library(metagMisc)
source("~/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/phyloseq_funcs.R", local = TRUE) ### load in functions for analysis

colorlist<-c("#1b6ca8", "#7CBBDE", "#438876", "#61C3AA", "#C06029", "#F6B38C", "#B88706", "#E8B120", "#CD4262", "#E18DA0")

## Reading in Artifacts (qza)
## HMP.feature.table<-read_qza("HMP.feature-table.qza")
## HMP.feature.tax.table<-read_qza("HMP.feature-tax-table.qza")

# physeq <- readRDS("/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/raw_HMP_FSJ_physeq.rds")

# ps.FSJ.cleaned <- readRDS("/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/noHMP.noSapien.rds")

# physeq.noncontam <- readRDS() # all samples, but contaminants removed

# physeq.noncontam.prev1 <- prune_taxa(!contamdf.prev1$contaminant, physeq)


# OR create your Phyloseq with new metadata
physeq<-qza_to_phyloseq(
  features= "../removing_HMP/HMP-table-no-hsapiens-exact.qza",
  taxonomy = "../removing_HMP/HMP-feature-tax-table-2022.qza",
  metadata = "../removing_HMP/FSJ_HMP_Metadata_table_2022.txt") #5736 unfiltered taxa

allSamples <- read.delim("../removing_HMP/FSJ_HMP_Metadata_table_2022.txt", header = T)
FSJsamples <- allSamples[allSamples$Sample_or_Control == "True Sample", ]
FSJsamplesID <- FSJsamples$sample_name

################## MAKE/LOAD IN YOUR MULTIPLE DIFFERENT PHYLOSEQ OBJECTS ########################
################## total OTU abundance is >= 5 reads it'll be preserved too
physeq <- phyloseq_filter_prevalence(physeq, prev.trh = 0, abund.trh = 5, threshold_condition = "AND", abund.type = "total")  # 4447

# Keep only the most abundant five phyla. (FOR OTU ORDINATION)
phylum.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
ps1 = prune_taxa((tax_table(physeq)[, "Phylum"] %in% top5phyla), physeq) #pruned to the top 5 phylum

physeq.noncontam # HMP + FSJ with contaminants removed
FSJphyseq = ps_filter(physeq, Sample_or_Control != "Control Sample", .keep_all_taxa = F) # filter out the HMP data -- but not CLEANED
FSJphyseq_cleaned # only FSJ, HMP filtered out
fracFSJ = transform_sample_counts(FSJphyseq_cleaned , function(x) x / sum(x)) #get proportion

FSJ.P1 # FSJ Prev1 strict filtering!
FSJ.P1 <- phyloseq_filter_prevalence(FSJ.P1, prev.trh = 0, abund.trh = 5, threshold_condition = "AND", abund.type = "total") 
FSJ.P1
################## check that you have no empty taxa (where no sample contains it -- rowSums) #3817 taxa, N=291
head(table(rowSums(ps1@otu_table))) #looks good
################## check that you have no samples without taxa (colSums)
head(table(colSums(ps1@otu_table))) # 33 individuals have no taxa (N=1904)
head(table(colSums(physeq@otu_table))) # 22 individuals
head(table(colSums(FSJ.P1@otu_table)))

# get sample names that don't have taxa
ps1 <- phyloseq_richness_filter(ps1, mintaxa = 1) # now 1871 samples
physeq <- phyloseq_richness_filter(physeq, mintaxa = 1) 
head(table(colSums(ps1@otu_table))) #should not have a zero column anymore.
head(table(colSums(physeq@otu_table))) 
head(table(colSums(FSJ.P1@otu_table))) 

##### Label HUMAN HMP samples
physeq@sam_data[["sample_type"]]<-as.factor(physeq@sam_data[["sample_type"]])
sample_data(physeq)$human <- get_variable(physeq, "sample_type") %in% c("HMP_oral", "HMP_nasal", "HMP_skin")
sample_data(physeq)$FSJ <- get_variable(physeq, "sample_type") %in% c("ABS_Contemporary", "ABS_Historic", "JDSP_Contemporary", "JDSP_Historic", "ONF_Contemporary", "ONF_Historic", "PLE_Contemporary", "PLE_Historic", "SPSP_Historic", "SPSP_Contemporary")

#reorder our populations to plot them 
sample_data(physeq)$sample_type <- factor(sample_data(physeq)$sample_type, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary", "HMP_nasal", "HMP_oral", "HMP_skin"))
sample_data(physeq) <- sample_data(physeq)[order(sample_data(physeq)$sample_type), ]


### FSJ ONLY
sample_data(FSJ.P1)$sample_type <- factor(sample_data(FSJ.P1)$sample_type, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
sample_data(FSJ.P1) <- sample_data(FSJ.P1)[order(sample_data(FSJ.P1)$sample_type), ]

############## Plot how many taxa found in each FSJ before Decontam ############

c.df <- as.data.frame(cbind(colSums(FSJphyseq@otu_table), FSJsamples$sample_name, FSJsamples$sample_type))
colnames(c.df) <- c("NumOTUs", "Sample", "Category")
c.df$NumOTUs <- as.numeric(c.df$NumOTUs)

ggplot(c.df, aes(x=reorder(Sample,-NumOTUs), y=NumOTUs, color=Category))+
  geom_point()+
  xlab("Samples")+ ylab("Number of OTUs found") +
  #scale_color_manual(values=c("#f07605", "#020fc9", "#54B435", "white"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme( # remove the vertical grid lines 
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.1, color="black" )) +
  theme(legend.title = element_text(size=14))+
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size = 14), axis.ticks.x=element_blank())+
  theme(axis.title.y = element_text(vjust=5))+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(plot.margin = margin(12, 12, 10, 22))

################## Get count of TAXA COUNT per sample

####make absence and presence OTU table
FSJphyseq@otu_table[FSJphyseq@otu_table == 0.00000]<-0 #absent
FSJphyseq@otu_table[FSJphyseq@otu_table != 0.00000]<- 1 #present

c.df <- as.data.frame(cbind(colSums(FSJphyseq@otu_table), FSJsamples$sample_name, FSJsamples$sample_type))
colnames(c.df) <- c("NumTaxa", "Sample", "Category")
c.df$NumTaxa <- as.numeric(c.df$NumTaxa)

#reorder our populations to plot them 
c.df$Category<- factor(c.df$Category, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
c.df <- c.df[order(c.df$Category), ]

ggplot(c.df, aes(x=Category, y=NumTaxa, fill=Category))+
  geom_boxplot()+
  ylab("Number of taxa found") +
  scale_fill_manual(values=c(colorlist))+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("") + ylab("Total number of ROHs")+
  theme(legend.title=element_blank(), legend.text = element_text(size=14), legend.background =element_blank()) +
  theme(legend.position = "")+
  theme(legend.text = element_text(size=18))+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18, vjust=1, hjust=0), axis.text.y=element_text(size = 18))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))+ #trbl
  #ggtitle("ROH lengths per individual in ABS and PLE")+
  scale_x_discrete(labels=c("ABS 1992", "ABS 2017", "ONF 1992", "ONF 2017", "JDSP 1992", "JDSP 2017", "PLE 1992", "PLE 2017", "SPSP 2004", "SPSP 2017"))

# Get summary by population
tapply(c.df$NumTaxa, c.df$Category, summary)


ggplot(c.df, aes(NumTaxa))+
  geom_histogram(fill="steelblue")+
  ylab("Frequency (count of samples)") + xlab("Number of distinct taxa found")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank(), legend.text = element_text(size=14), legend.background =element_blank()) +
  theme(legend.position = "")+
  theme(legend.text = element_text(size=18))+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18, vjust=1, hjust=0), axis.text.y=element_text(size = 18))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))




##### CONTAMINANTS FILTERED OUT
####make absence and presence OTU table
physeq.noncontam@otu_table[physeq.noncontam@otu_table == 0.00000]<-0 #absent
physeq.noncontam@otu_table[physeq.noncontam@otu_table != 0.00000]<- 1 #present

c.df <- as.data.frame(cbind(colSums(physeq.noncontam@otu_table), FSJsamples$sample_name, FSJsamples$sample_type))
colnames(c.df) <- c("NumTaxa", "Sample", "Category")
c.df$NumTaxa <- as.numeric(c.df$NumTaxa)

#reorder our populations to plot them 
c.df$Category<- factor(c.df$Category, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
c.df <- c.df[order(c.df$Category), ]

ggplot(c.df, aes(x=Category, y=NumTaxa, fill=Category))+
  geom_boxplot()+
  ylab("Number of taxa found") +
  scale_fill_manual(values=c(colorlist))+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("") + ylab("Total number of ROHs")+
  theme(legend.title=element_blank(), legend.text = element_text(size=14), legend.background =element_blank()) +
  theme(legend.position = "")+
  theme(legend.text = element_text(size=18))+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18, vjust=1, hjust=0), axis.text.y=element_text(size = 18))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))+ #trbl
  #ggtitle("ROH lengths per individual in ABS and PLE")+
  scale_x_discrete(labels=c("ABS 1992", "ABS 2017", "ONF 1992", "ONF 2017", "JDSP 1992", "JDSP 2017", "PLE 1992", "PLE 2017", "SPSP 2004", "SPSP 2017"))

# Get summary by population
tapply(c.df$NumTaxa, c.df$Category, summary)



############# Get a bar plot of the proportions of the 5 top phyla, and the remainder as "others"

y1 <- tax_glom(FSJphyseq_cleaned, taxrank = 'Phylum') # agglomerate taxa
(y2 = merge_samples(y1, "sample_type")) # merge samples on sample variable of interest
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Phylum[y4$Abundance < 0.01] <- "Phylum < 1% abund." #rename genera with < 1% abundance

#set color palette to accommodate the number of genera
colourCount = length(unique(y4$Phylum))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#plot
y4$Sample<- factor(y4$Sample, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
y4 <- y4[order(y4$Sample), ]

unique(y4$Phylum)
y4$Phylum<- factor(y4$Phylum, levels = c("Proteobacteria","Actinobacteria","Firmicutes","Cyanobacteria","Euryarchaeota","Bacteroidetes","Phylum < 1% abund."))
y4 <- y4[order(y4$Phylum), ]

p <- ggplot(data=y4, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=getPalette(colourCount)) + guides(fill=guide_legend(nrow=length(unique(y4$Phylum)))) +theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 15, vjust=0.5, size=18), axis.text.y=element_text(size = 18))

p
ggsave(p, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/phylum_stacked_noContam-FSJ.pdf", width=15, height = 8, units = "in", device="pdf", dpi=800)


## PLOT GENUS
y1 <- tax_glom(FSJphyseq_cleaned, taxrank = 'Genus') # agglomerate taxa
(y2 = merge_samples(y1, "sample_type")) # merge samples on sample variable of interest
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Genus <- as.character(y4$Genus) #convert to character
y4$Genus[y4$Abundance < 0.01] <- "Genera < 1% abund." #rename genera with < 1% abundance

#set color palette to accommodate the number of genera
colourCount = length(unique(y4$Genus))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#plot
y4$Sample<- factor(y4$Sample, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
y4 <- y4[order(y4$Sample), ]

y4$Genus<- factor(y4$Genus, levels = c("Alcanivorax","Streptomyces","Paraburkholderia","Nocardioides","Devosia","Achromobacter","Rathayibacter","Xanthomonas","Luteimonas","Frankia","Haloarcula","Micromonospora","Roseococcus","Protaetiibacter","Pseudonocardia","Methylomonas","Pseudomonas","Leptolyngbya", "Ornithinimicrobium","Brevundimonas","Serratia","Stenotrophomonas","Microbacterium","Bacillus","Campylobacter", "Genera < 1% abund."))
y4 <- y4[order(y4$Genus), ]

p <- ggplot(data=y4, aes(x=Sample, y=Abundance, fill=Genus)) + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + guides(fill=guide_legend(nrow=length(unique(y4$Genus)))) +theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 15, vjust=0.5, size=18), axis.text.y=element_text(size = 18))
p
ggsave(p, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/Genera_stacked_noContam-FSJ.pdf", width=15, height = 8, units = "in", device="pdf", dpi=800)


######################## PLOT FAMILY ##########################
y1 <- tax_glom(FSJphyseq_cleaned, taxrank = 'Family') # agglomerate taxa
(y2 = merge_samples(y1, "sample_type")) # merge samples on sample variable of interest
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Family <- as.character(y4$Family) #convert to character
y4$Family[y4$Abundance < 0.01] <- "Families < 1% abund." #rename genera with < 1% abundance

#set color palette to accommodate the number of genera
colourCount = length(unique(y4$Family))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#plot
y4$Sample<- factor(y4$Sample, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
y4 <- y4[order(y4$Sample), ]


# get family vector
vec<- as.character(unique(y4$Family))
paste(vec, sep=",", collapse = ",")


"Alcanivoracaceae", "Streptomycetaceae","Burkholderiaceae","Microbacteriaceae","Xanthomonadaceae" ,"Nocardioidaceae"       "Hyphomicrobiaceae","Alcaligenaceae","Pseudonocardiaceae","Acetobacteraceae",
"Haloarculaceae","Frankiaceae","Micromonosporaceae","Methylococcaceae","Pseudomonadaceae" ,   
"Yersiniaceae","Rhodobacteraceae","Leptolyngbyaceae","Ornithinimicrobiaceae",
"Caulobacteraceae","Sphingomonadaceae","Bacillaceae","Campylobacteraceae", "Families < 1% abund."

y4$Family<- factor(y4$Family, levels = c("Alcanivoracaceae","Streptomycetaceae","Burkholderiaceae","Microbacteriaceae","Xanthomonadaceae" ,"Nocardioidaceae","Hyphomicrobiaceae","Alcaligenaceae","Pseudonocardiaceae","Acetobacteraceae","Haloarculaceae","Frankiaceae","Micromonosporaceae","Methylococcaceae","Pseudomonadaceae","Yersiniaceae","Rhodobacteraceae","Leptolyngbyaceae","Ornithinimicrobiaceae","Caulobacteraceae","Sphingomonadaceae","Bacillaceae","Campylobacteraceae", "Families < 1% abund."))

y4 <- y4[order(y4$Family), ]

p <- ggplot(data=y4, aes(x=Sample, y=Abundance, fill=Family)) + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCount)) + guides(fill=guide_legend(nrow=length(unique(y4$Family)))) +theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 15, vjust=0.5, size=18), axis.text.y=element_text(size = 18))
p
ggsave(p, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/Family_stacked_noContam-FSJ.pdf", width=15, height = 8, units = "in", device="pdf", dpi=800)



#plot_bar(ps1, fill="Phylum") ## this take very long, be careful. dont do with Genus or on.
#plot_bar(ps1, x="sample_type", fill="Family")


########## ALPHA BETA DIVERSITY ##############

# prune OTUs that are not present in any of the samples
GP <- prune_taxa(taxa_sums(FSJphyseq_cleaned) > 0, FSJphyseq_cleaned)

# calculate several metrics of diversity
alpha <- estimate_richness(GP)

## add meta data to your alpha diversity
meta <- read.delim("../removing_HMP/FSJ_Metadata.txt", header = T, as.is=T, sep = "\t")

AMeta <- as.data.frame(cbind(alpha, meta))

# reorder
AMeta$sample_type<- factor(AMeta$sample_type, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
AMeta <- AMeta[order(AMeta$sample_type), ]

# get all pairwise comparisons in GGPLOT2 -- compare means
comp <- list(c("ONF_historic", "ONF_contemporary"), c("ABS_historic", "ABS_contemporary"), c("PLE_historic", "PLE_contemporary"), c("JDSP_historic", "JDSP_contemporary"), c("SPSP_historic", "SPSP_contemporary")) # within populations, across years

pwc <- AMeta %>%
  pairwise_wilcox_test(Shannon ~ sample_type)
pwc
pwc <- pwc %>% add_xy_position(x = "sample_type")


a <- ggboxplot(AMeta, x="sample_type", y="Shannon", fill="sample_type", palette = c(colorlist))+
  theme_minimal()+ ylab("Shannon Index")+
  stat_pvalue_manual(pwc, hide.ns = TRUE, label = "p.adj.signif") +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )+
  theme(legend.position = "none")+
  geom_jitter(size=1.5, shape=19, position=position_jitter(0))+
  scale_fill_manual(values = colorlist)+ theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 15, vjust=0.5, size=18), axis.text.y=element_text(size = 18))
a
a +  theme(legend.position = "none")

### Plot Simpsons
a <- ggplot(AMeta, aes(x=sample_type, y=Simpson, fill=sample_type))+
  geom_boxplot()+
  theme_minimal()+ ylab("Simpson Index")+
  scale_fill_manual(values = colorlist)+ theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  geom_jitter(position=position_jitter(0))+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 15, vjust=0.5, size=18), axis.text.y=element_text(size = 18))+
  theme(plot.margin = margin(15, 22, 15, 22))

a +  theme(legend.position = "none")

############################### HET VERSUS ALPHA ##################################
all_het <- read.delim("/Users/tramnguyen/Box Sync/PhD/Projects/FSJ-WGS20x/2022-July/HET/breeders_snps_hwe_noZ_genomind15_noLD.het", as.is=T, header=T, sep="")
all_het <- all_het %>% mutate(Obs.Het = 1-(O.HOM./N.NM.))
all_het$IID <- paste0(all_het$IID, "_stdb") # change sample names to match kraken samples
colnames(all_het)[2] <- "sample_name"

# merge het and alpha diversity
aHet <- merge(AMeta, all_het, by="sample_name")

# plot shannon versus het
library("ggpubr")
ggscatter(aHet, x = "Obs.Het", y = "Shannon", color = "sample_type", 
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")

ONF <- subset(aHet, aHet$sample_type == "ONF_Historic" | aHet$sample_type == "ONF_Contemporary",)

ggscatter(ONF, x = "Obs.Het", y = "Shannon",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")

PLE <- subset(aHet, aHet$sample_type == "PLE_Historic" | aHet$sample_type == "PLE_Contemporary",)

ggscatter(PLE, x = "Obs.Het", y = "Shannon",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


JDSP <- subset(aHet, aHet$sample_type == "JDSP_Historic" | aHet$sample_type == "JDSP_Contemporary",)

ggscatter(JDSP, x = "Obs.Het", y = "Shannon",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


ABS <- subset(aHet, aHet$sample_type == "ABS_Historic" | aHet$sample_type == "ABS_Contemporary",)

ggscatter(ABS, x = "Obs.Het", y = "Shannon",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


SPSP <- subset(aHet, aHet$sample_type == "SPSP_Historic" | aHet$sample_type == "SPSP_Contemporary",)

ggscatter(SPSP, x = "Obs.Het", y = "Shannon",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


cor.test(aHet$Obs.Het, aHet$Shannon, method = c("pearson"))


plot_richness(GP, x="sample_type", measures=c("Shannon", "Simpson"))



################### SIMPSON #####################

library("ggpubr")
ggscatter(aHet, x = "Obs.Het", y = "Simpson", color = "sample_type", 
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")

ONF <- subset(aHet, aHet$sample_type == "ONF_Historic" | aHet$sample_type == "ONF_Contemporary",)

ggscatter(ONF, x = "Obs.Het", y = "Simpson",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")

PLE <- subset(aHet, aHet$sample_type == "PLE_Historic" | aHet$sample_type == "PLE_Contemporary",)

ggscatter(PLE, x = "Obs.Het", y = "Simpson",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


JDSP <- subset(aHet, aHet$sample_type == "JDSP_Historic" | aHet$sample_type == "JDSP_Contemporary",)

ggscatter(JDSP, x = "Obs.Het", y = "Simpson",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


ABS <- subset(aHet, aHet$sample_type == "ABS_Historic" | aHet$sample_type == "ABS_Contemporary",)

ggscatter(ABS, x = "Obs.Het", y = "Simpson",
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")


cor.test(aHet$Obs.Het, aHet$Shannon, method = c("pearson"))


r <- ggplot(aes(x = Obs.Het, y = Shannon, color = sample_type, group = 1), data = aHet) +
  geom_point(size=4, alpha=0.8) + scale_color_manual(values=c(colorlist)) + 
  geom_smooth(color="black", method="lm", se=T, linetype="solid", fullrange=TRUE, ) +
  theme_minimal()+ ylab("Shannon Index")+ xlab("Observed host heterozygosity")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18), axis.text.y=element_text(size = 18))+
  ggtitle("Het versus Shannon using FSJ_Cleaned")
r
ggsave(r, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/HET_shannon1.pdf", width=9, height = 6, units = "in", device="pdf", dpi=800)




r <- ggplot(aes(x = Obs.Het, y = Simpson, color = sample_type, group = 1), data = aHet) +
  geom_point(size=4, alpha=0.8) + scale_color_manual(values=c(colorlist)) + 
  geom_smooth(color="black", method="lm", se=T, linetype="solid", fullrange=TRUE, ) +
  theme_minimal()+ ylab("Simpson Index")+ xlab("Observed host heterozygosity")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18), axis.text.y=element_text(size = 18))+
  ggtitle("Het versus Simpson using FSJ_Cleaned")
r
ggsave(r, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/HET_simpson1.pdf", width=9, height = 6, units = "in", device="pdf", dpi=800)

cor.test(aHet$Obs.Het, aHet$Simpson, method = c("pearson"))


############################################ FROH VERSUS ALPHA #########################################

all_segments<-read.delim(file = "/Users/tramnguyen/Box Sync/PhD/Projects/FSJ-WGS20x/2022-July/ROH/N268_noZ_hwe_genomind15.hom", header = T,  as.is = T, sep = "")
roh_inds<-read.delim(file = "/Users/tramnguyen/Box Sync/PhD/Projects/FSJ-WGS20x/2022-July/ROH/N268_noZ_hwe_genomind15.hom.indiv", header = T,  as.is = T, sep = "")
#pops<-read.delim("/Users/tramnguyen/Box Sync/PhD/Projects/FSJ-WGS20x/Scripts/pop_list_analyses.txt", header = F)
contig_lengths<-read.delim("~/Box Sync/PhD/Projects/FSJ-WGS20x/Updated2021/Filter7/Variant_QC/scaffold_lengths.txt", header = F) #Z = ScYP8k3_10 = 75605511bp
# sum(contig_lengths$V2) = 1060970510
# contigs minus Z = 985364999 (in KB = 985364.999)


roh_inds_prop <- roh_inds %>% mutate(FROH = KB/985364999)


roh_inds_prop$IID <- paste0(roh_inds_prop$IID, "_stdb") # change sample names to match kraken samples
colnames(roh_inds_prop)[2] <- "sample_name"
head(roh_inds_prop)

# merge het and alpha diversity
aFROH <- merge(AMeta, roh_inds_prop, by="sample_name")



##########  PLOT ###########

library(tidyverse)
library(ggpmisc)
my.formula = y ~ x

r <- ggplot(aes(x = FROH, y = Shannon, color = sample_type, group = 1), data = aFROH) +
  geom_point(size=4, alpha=0.8) + scale_color_manual(values=c(colorlist)) + 
  geom_smooth(color="black", method="lm", se=T, linetype="solid", fullrange=TRUE, ) +
  theme_minimal()+ ylab("Shannon Index")+ xlab("Proportion of genome in ROH (Inbreeding)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18), axis.text.y=element_text(size = 18))+
  ggtitle("FROH versus Shannon using FSJ_Cleaned")
r
ggsave(r, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/FROH_shannon1.pdf", width=9, height = 6, units = "in", device="pdf", dpi=800)


ggscatter(aFROH, x = "FROH", y = "Shannon", color = "sample_type", 
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black", fullrange=T)

cor.test(aFROH$FROH, aFROH$Shannon, method = c("pearson"))

###### SIMPSON #####

r <- ggplot(aes(x = FROH, y = Simpson, color = sample_type, group = 1), data = aFROH) +
  geom_point(size=4, alpha=0.8) + scale_color_manual(values=c(colorlist)) + 
  geom_smooth(color="black", method="lm", se=T, linetype="solid", fullrange=TRUE, ) +
  theme_minimal()+ ylab("Simpson Index")+ xlab("Proportion of genome in ROH (Inbreeding)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18), axis.text.y=element_text(size = 18))+
  ggtitle("FROH versus Shannon using FSJ_Cleaned")
r
ggsave(r, filename = "~/Box Sync/PhD/Projects/Kraken_Pathogen/phyloseq/FROH_simpson1.pdf", width=9, height = 6, units = "in", device="pdf", dpi=800)

cor.test(aFROH$FROH, aFROH$Simpson, method = c("pearson"))



########################## ALPHA DIVERSITY WITH JUST JDSP ADULT SAMPLES COLLECTED BY ME ####################
#### get a list of JDSP adults that I banded in the same sequencing bath
mylist <- c("JDSP027_stdb", "JDSP089_stdb", "JDSP090_stdb", "JDSP130_stdb", "JDSP138_stdb", "JDSP139_stdb", "JDSP140_stdb", "JDSP141_stdb", "JDSP142_stdb", "JDSP143_stdb", "JDSP144_stdb", "JDSP145_stdb", "JDSP146_stdb", "JDSP147_stdb", "JDSP148_stdb", "JDSP149_stdb", "JDSP150_stdb", "JDSP151_stdb", "JDSP902_stdb", "JDSP909_stdb")

mylist <- c("JDSP027_stdb", "JDSP089_stdb", "JDSP090_stdb", "JDSP130_stdb", "JDSP138_stdb", "JDSP139_stdb", "JDSP140_stdb", "JDSP141_stdb", "JDSP142_stdb", "JDSP143_stdb", "JDSP144_stdb", "JDSP146_stdb", "JDSP147_stdb", "JDSP149_stdb", "JDSP150_stdb", "JDSP151_stdb", "JDSP909_stdb") #remove the three that separated out in PCA space for now


jdsp_TN<-prune_samples(mylist,FSJphyseq_cleaned)
jdsp_TN <- prune_taxa(taxa_sums(jdsp_TN)>0, jdsp_TN)

sample_data(jdsp_TN)$sample_name<-mylist

jdsp_TN.ord <- ordinate(jdsp_TN)
plot_ordination(jdsp_TN, jdsp_TN.ord, type="samples", color="sample_type", title="taxa", label="sample_name")+ geom_label_repel(aes(label = sample_name), size = 3)


all_het <- read.delim("/Users/tramnguyen/Box Sync/PhD/Projects/FSJ-WGS20x/2022-July/HET/breeders_snps_hwe_noZ_genomind15_noLD.het", as.is=T, header=T, sep="")
all_het <- all_het %>% mutate(Obs.Het = 1-(O.HOM./N.NM.))
all_het$IID <- paste0(all_het$IID, "_stdb") # change sample names to match kraken samples
colnames(all_het)[2] <- "sample_name"

# merge het and alpha diversity
aHet <- merge(AMeta, all_het, by="sample_name")

c.Het<-subset(aHet, aHet$sample_name %in% mylist)

# plot shannon versus het
library("ggpubr")
ggscatter(c.Het, x = "Obs.Het", y = "Shannon", color = "sample_type", 
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black")




############# What about inbreeding ###############

c.FROH<-subset(aFROH, aFROH$sample_name %in% mylist)

library(tidyverse)
library(ggpmisc)
my.formula = y ~ x

r <- ggplot(aes(x = FROH, y = Shannon, color = sample_type, group = 1), data = c.FROH) +
  geom_point(size=4, alpha=0.8) + scale_color_manual(values=c(colorlist)) + 
  geom_smooth(color="black", method="lm", se=T, linetype="solid", fullrange=TRUE, ) +
  theme_minimal()+ ylab("Shannon Index")+ xlab("Proportion of genome in ROH (Inbreeding)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust=0.5, size=24))+
  theme(legend.title=element_text(size=18), legend.text = element_text(size=16), legend.background =element_blank()) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=18), axis.text.y=element_text(size = 18))+
  ggtitle("FROH versus Shannon using FSJ_Cleaned")
r


ggscatter(aFROH, x = "FROH", y = "Shannon", color = "sample_type", 
          add = "reg.line", conf.int = T, size = 2,
          cor.coef = TRUE, cor.method = "pearson", palette = c(colorlist), xlab = "", ylab="")+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth", color = "black", fill="black", fullrange=T)

cor.test(aFROH$FROH, aFROH$Shannon, method = c("pearson"))


##################### Beta diversity #######################

######################### COMPARE PAIRWISE IBD WITH OTU DISSIMILARITY #########################
# https://rpubs.com/lconteville/714853
pIBD <- read.delim("/Users/tramnguyen/Box Sync/PhD/Projects/FSJ-WGS20x/2022-July/IBD/N268_noZ_hwe_genomind15_noLD.genome", header = T,  as.is = T, sep = "")

pIBD$IID1 <- paste0(pIBD$IID1, "_stdb")
pIBD$IID2 <- paste0(pIBD$IID2, "_stdb")

pIBD_samples<-sort(unique(c(pIBD$IID1, pIBD$IID2))) # get vector of samples in kraken format that have genetic data for


### Make a matrix of your pairwise IBD
library(reshape2)

## Informs acast that you have all these different values and it should make it the levels
sample.list <- with(pIBD, sort(unique(c(IID1, IID2))))
pIBD$IID1 <- with(pIBD, factor(IID1, levels = sample.list))
pIBD$IID2 <- with(pIBD, factor(IID2, levels = sample.list))

# check that there are no NAs
#range(df1$PI_HAT)
#0.0000 0.7327

#convert long-to-wide
x <- acast(pIBD, IID1 ~ IID2, value.var = "PI_HAT", fill=0, drop=F) ## filling = 1 because missing values are the within same sample IBD that is missing
# fill = NA set everything else that is missing to NA
x  = x+t(x) ### if you had this originally as NA, it couldn't add them
diag(x) <- 1 ## set diagonals to 1
x[1:5,1:5]
isSymmetric(x)

### PHYLOSEQ MANIPULATION 
# extract only the 268 individuals that you have popgen data for (pairwise IBD)

GenoFSJ= prune_samples(pIBD_samples, FSJphyseq_cleaned)
fracFSJ = transform_sample_counts(GenoFSJ, function(x) x / sum(x) * 100) 
head(otu_table(fracFSJ)[,1:6])

## calculate dissimilarity between samples
fsj_bray <- phyloseq::distance(fracFSJ, method = "bray")
fsj_bray <- as.matrix(fsj_bray)
head(fsj_bray)[,1:6]


mylong <- data.frame(IBD=x[upper.tri(x)], bray=fsj_bray[upper.tri(fsj_bray)], stringsAsFactors = F)

plot(mylong$IBD, mylong$bray)

## do a Mantel test
dim(x)
dim(fsj_bray)
mantel(x,fsj_bray,method = "pearson") 
#Mantel statistic r: -0.02949 Significance: 0.991


## The “fsj_bray” matrix presents the distance between all samples, but since we wanna generate a boxplot with distances considering the metagenomes in each group separately, we need to filter this matrix. With the code below we end with a dataframe “df.bray” storing Bray Curtis Distances between metagenomes from the same groups.

sub_dist <- list()
groups_all <- sample_data(fracFSJ)$sample_type

for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- sample_names(fracFSJ)[row_group]
  sub_dist[[group]] <- fsj_bray[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))

head(df.bray)


ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(alpha=0.6) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))

# To generate the PcoA plot, first get the ordination result with the “ordinate” command from phyloseq. It requires a phyloseq object, and it accepts diverse methods and distances. Next, generate the plot using the “plot_ordination()” function.
ord = ordinate(fracFSJ, method="PCoA", distance = "bray")

plot_ordination(fracFSJ, ord, color = "sample_type", shape="sample_type") + 
  geom_point(size=4) + 
  stat_ellipse(aes(group=sample_type))

#To test if the Control and the Treatment Groups are significantly different from each other we can run a Permanova test using the adonis function from the vegan package.

samples <- data.frame(sample_data(fracFSJ))
adonis(fsj_bray ~ sample_type, data = samples)

#This result means that the adonis test is significant, therefore the null hypothesis in which the Control/Treatment have the same centroid was reject.



### Scatterplot of pairwise IBD and bray curtis 
### each point will be one pairwise comparison


# Calculate distances
DistBC = distance(FSJphyseq_cleaned, method = "bray")
dist_bc <- as.matrix(vegdist(FSJphyseq_cleaned@otu_table, method = "bray")) 
dist_bc[1:5, 1:5]

# calculate PCOA using Phyloseq package
pcoa_bc = ordinate(FSJphyseq_cleaned@otu_table, "PCoA", "bray") 
pcoa_bc

plot_ordination(FSJphyseq_cleaned, pcoa_bc, color = "Age") + 
  geom_point(size = 3) 




# calculate Bray-Curtis distance using the vegan package