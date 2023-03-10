## Plot PCAs Ordination 
## Dec 2022

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

################## MAKE YOUR FOUR DIFFERENT PHYLOSEQ OBJECTS ########################
################## total OTU abundance is >= 5 reads it'll be preserved too
physeq <- phyloseq_filter_prevalence(physeq, prev.trh = 0, abund.trh = 5, threshold_condition = "AND", abund.type = "total")  # 4447

# Keep only the most abundant five phyla. (FOR OTU ORDINATION)
phylum.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
ps1 = prune_taxa((tax_table(physeq)[, "Phylum"] %in% top5phyla), physeq) #pruned to the top 5 phylum

physeq.noncontam # HMP + FSJ with contaminants removed
FSJphyseq = ps_filter(physeq, Sample_or_Control != "Control Sample", .keep_all_taxa = F) # filter out the HMP data
FSJphyseq
FSJphyseq_cleaned # only FSJ, HMP filtered out


FSJ.P1 # FSJ Prev1 strict filtering!
FSJ.P1 <- phyloseq_filter_prevalence(FSJ.P1, prev.trh = 0, abund.trh = 5, threshold_condition = "AND", abund.type = "total") 
FSJ.P1

################## check that you have no empty taxa (where no sample contains it -- rowSums) #3817 taxa, N=291
head(table(rowSums(FSJphyseq@otu_table))) #looks good
################## check that you have no samples without taxa (colSums)
head(table(colSums(FSJphyseq@otu_table))) # 33 individuals have no taxa (N=1904)
head(table(colSums(physeq@otu_table))) # 22 individuals
head(table(colSums(physeq.noncontam@otu_table)))
head(table(colSums(FSJphyseq_cleaned@otu_table)))


# get sample names that don't have taxa
ps1 <- phyloseq_richness_filter(ps1, mintaxa = 1) # now 1871 samples
physeq <- phyloseq_richness_filter(physeq, mintaxa = 1) 
physeq.noncontam <- phyloseq_richness_filter(physeq.noncontam, mintaxa = 1) 
head(table(colSums(ps1@otu_table))) #should not have a zero column anymore.
head(table(colSums(physeq@otu_table))) 
head(table(colSums(physeq.noncontam@otu_table))) 

##### Label HUMAN HMP samples
physeq.noncontam@sam_data[["sample_type"]]<-as.factor(physeq.noncontam@sam_data[["sample_type"]])
sample_data(physeq.noncontam)$human <- get_variable(physeq.noncontam, "sample_type") %in% c("HMP_oral", "HMP_nasal", "HMP_skin")
sample_data(physeq.noncontam)$FSJ <- get_variable(physeq.noncontam, "sample_type") %in% c("ABS_Contemporary", "ABS_Historic", "JDSP_Contemporary", "JDSP_Historic", "ONF_Contemporary", "ONF_Historic", "PLE_Contemporary", "PLE_Historic", "SPSP_Historic", "SPSP_Contemporary")

#reorder our populations to plot them 
sample_data(physeq.noncontam)$sample_type <- factor(sample_data(physeq.noncontam)$sample_type, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary", "HMP_nasal", "HMP_oral", "HMP_skin"))
sample_data(physeq.noncontam) <- sample_data(physeq.noncontam)[order(sample_data(physeq.noncontam)$sample_type), ]


sample_data(FSJphyseq_cleaned)$sample_type <- factor(sample_data(FSJphyseq_cleaned)$sample_type, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
sample_data(FSJphyseq_cleaned) <- sample_data(FSJphyseq_cleaned)[order(sample_data(FSJphyseq_cleaned)$sample_type), ]


########################### ORDINATION PLOTS WITHOUT DECONTAM ##########################
################## Four main ordination plots #####################
# 1. Just OTUs

GP.ord <- ordinate(ps1, "NMDS", "bray", "PCoA")
p1 = plot_ordination(ps1, GP.ord, type="taxa", color="Phylum", title="PCoA of the five most abundant Phyla")
p1 + geom_point(size=2, alpha=0.9) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray 93"),
        panel.grid.minor = element_line(color="gray 94"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x = element_text(size=14), axis.text.y=element_text(size = 14))+
  theme(plot.title=element_text(hjust=0.5, size=18))
print(p1)

#p1 + facet_wrap(~Phylum, 3)



# 2. Just samples

samp.ord <- ordinate(FSJphyseq_cleaned,"MDS", "bray", "PCoA")
p1 = plot_ordination(FSJphyseq_cleaned, samp.ord, type="samples", color="sample_type", title="PCoA of FSJ samples after contaminant filtering")+ geom_point(size=3, alpha=0.9) +
  scale_color_manual(values=c(colorlist, "#371B58", "#5B4B8A", "#916ec4"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray 93"),
        panel.grid.minor = element_line(color="gray 94"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust = 0.6), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=5))+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(plot.margin = margin(12, 12, 10, 22))+
  theme(plot.title=element_text(hjust=0.5, vjust=2.5, size=18))

p1

ggsave(p1, file="/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/PCoA_FSJ_cleaned.pdf", width = 12, height=8, units = "in", device="pdf", dpi=800) # with physeq.noncontam objects

ggsave(p1, file="/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/PCoA_all_samp_cleaned.pdf", width = 12, height=8, units = "in", device="pdf", dpi=800) # with physeq.noncontam objects
ggsave(p1, file="/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/PCoA_all_samp_NoDecontamFilter.pdf", width = 12, height=8, units = "in", device="pdf", dpi=800) # with physeq object


################## COLORED BY COLLECTOR
samp.ord <- ordinate(FSJphyseq_cleaned,"MDS", "bray", "PCoA")
p1 = plot_ordination(FSJphyseq_cleaned, samp.ord, type="samples", color="Collector", shape="Historic", title="PCoA of FSJ samples after contaminant filtering")+ geom_point(size=3, alpha=0.9) +
  scale_color_manual(values=c("#043b8c", "#c7511e", "#13752d"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray 93"),
        panel.grid.minor = element_line(color="gray 94"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust = 0.6), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=5))+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(plot.margin = margin(12, 12, 10, 22))+
  theme(plot.title=element_text(hjust=0.5, vjust=2.5, size=18))

p1

ggsave(p1, file="/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/PCoA_FSJ_cleaned_batch_Collector.pdf", width = 12, height=8, units = "in", device="pdf", dpi=800)


# 3. biplot graphic
# The plot_ordination function can also automatically create two different graphic layouts in which both the samples and OTUs are plotted together in one “biplot”. Note that this requires methods that are not intrinsically samples-only ordinations. For example, this doesn’t work with UniFrac/PCoA.

p3 = plot_ordination(ps1, GP.ord, type="biplot", color="sample_type", shape="Phylum", title="biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names = get_taxa_unique(physeq, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)


# 4. split graphic
# In this case the type="split" option can be helpful, in which the samples/OTUs are separated on two side-by-side panels…
p4 = plot_ordination(ps1, GP.ord, type="split", color="Phylum", label="sample_type", title="split") 
p4

gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
color.names <- levels(p4$data$Phylum)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols["samples"] <- "black"
p4 + scale_color_manual(values=p4cols)




########### add column for historic/contemporary
FSJphyseq_cleaned@sam_data[["sample_type"]]<-as.factor(physeq@sam_data[["sample_type"]])
sample_data(FSJphyseq_cleaned)$Historic <- get_variable(FSJphyseq_cleaned, "sample_type") %in% c("ABS_Historic", "ONF_Historic", "JDSP_Historic", "PLE_Historic", "SPSP_Historic")
sample_data(FSJphyseq_cleaned)$Contemporary <- get_variable(FSJphyseq_cleaned, "sample_type") %in% c("ABS_Contemporary", "JDSP_Contemporary", "ONF_Contemporary", "PLE_Contemporary", "SPSP_Contemporary")

#reorder our populations to plot them 
sample_data(FSJphyseq_cleaned)$sample_type <- factor(sample_data(FSJphyseq_cleaned)$sample_type, levels = c("ABS_Historic", "ABS_Contemporary", "ONF_Historic", "ONF_Contemporary", "JDSP_Historic", "JDSP_Contemporary", "PLE_Historic", "PLE_Contemporary", "SPSP_Historic", "SPSP_Contemporary"))
sample_data(FSJphyseq_cleaned) <- sample_data(FSJphyseq_cleaned)[order(sample_data(FSJphyseq_cleaned)$sample_type), ]

samp.ord <- ordinate(FSJphyseq_cleaned,"MDS", "bray", "PCoA")
p1 = plot_ordination(FSJphyseq_cleaned, samp.ord, type="samples", color="Collector", shape="Historic", title="PCoA of FSJ samples after contaminant filtering")+ geom_point(size=3, alpha=0.9) +
  scale_color_manual(values=c(colorlist, "#371B58", "#5B4B8A", "#916ec4"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray 93"),
        panel.grid.minor = element_line(color="gray 94"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust = 0.6), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=5))+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(plot.margin = margin(12, 12, 10, 22))+
  theme(plot.title=element_text(hjust=0.5, vjust=2.5, size=18))
p1

ggsave(p1, file="/Users/tramnguyen/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/PCoA_FSJ_cleaned_batch_collector.pdf", width = 12, height=8, units = "in", device="pdf", dpi=800)




# Optional pre-processing:
# Transform to even sampling depth. !!! We already have equal sampling depth...
# data is first transformed to relative abundance, creating the new GPr object, which is then filtered such that only OTUs with a mean greater than 10^-5 are kept
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

c.ps = transform_sample_counts(physeq, function(OTU) OTU/sum(OTU) )

# transforms abundance counts to fractional abundance
c.ps = transform_sample_counts(newPhyloObject, function(OTU) OTU/sum(OTU) )

# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
GP = filter_taxa(newPhyloObject, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Standardize abundances to the median sequencing depth
total = median(sample_sums(newPhyloObject))
standf = function(x, t=total) round(t * (x / sum(x)))
gps = transform_sample_counts(newPhyloObject, standf)

# Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)







 
 
 

## Command to filter out contaminants
 physeq.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
 physeq.noncontam

# clean zeroes
physeq.nozero <- prune_taxa(taxa_sums(physeq.noncontam)>0, physeq.noncontam)
physeq.nozero@sam_data[["sample_type"]]<-as.factor(physeq.nozero@sam_data[["sample_type"]])

## PCA of all FSJ and HMP -- Before Decontam pruning
physeq.all.nonzero <- prune_taxa(taxa_sums(physeq)>0, physeq)
physeq.all.nonzero@sam_data[["sample_type"]]<-as.factor(physeq.all.nonzero@sam_data[["sample_type"]])
sampleData(physeq.all.nonzero)$human <- getVariable(physeq.all.nonzero, "sample_type") %in% c("HMP_oral", "HMP_nasal", "HMP_skin")
sampleData(physeq.all.nonzero)$FSJ <- getVariable(physeq.all.nonzero, "sample_type") %in% c("ABS Contemporary", "ABS Historic", "JDSP Contemporary", "JDSP Historic", "ONF Contemporary", "ONF Historic", "PLE Contemporary", "PLE Historic", "SPSP Historic", "SPSP Contemporary")
physeq.all.ord <- ordinate(physeq.all.nonzero,"MDS")
p = plot_ordination(physeq.all.nonzero, physeq.all.ord, type="samples", color="sample_type", title="All samples PCA",shape="human")
p


#Creating Human and FSJ Sample Types
sampleData(physeq.nozero)$human <- getVariable(physeq.nozero, "sample_type") %in% c("HMP_oral", "HMP_nasal", "HMP_skin")
sampleData(physeq.nozero)$FSJ <- getVariable(physeq.nozero, "sample_type") %in% c("ABS Contemporary", "ABS Historic", "JDSP Contemporary", "JDSP Historic", "ONF Contemporary", "ONF Historic", "PLE Contemporary", "PLE Historic", "Savannas")

physeq.noncontam.ord <- ordinate(physeq.nozero,"MDS", "bray", "PCoA")
p = plot_ordination(physeq.nozero, physeq.noncontam.ord, type="samples", color="sample_type", title="All samples PCA",shape="human")
p<-p + theme(
  plot.title = element_text(color="black", size=20, hjust = 0.5),
  axis.title.x = element_text(color="black", size=14),
  axis.text.x = element_text(color="black", size=14),
  axis.text.y = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14),
  legend.text = element_text(size=12))
p
ggsave(p, file="All samples PCA.png", width = 9, height=6.5, units = "in", device="png")


library(ggrepel)

sample_data(FSJphyseq)$sample_name <- FSJsamplesID

physeqFSJ.ord <- ordinate(FSJphyseq,"MDS", "bray", "PCoA")
p1= plot_ordination(FSJphyseq, physeqFSJ.ord, type="samples", color="sample_type", title="All 291 FSJ samples PCA", label="label")
p1
p1 + geom_label_repel(aes(label = label), size = 3)


### label only ADULTS CAUGHT IN JDSP BY ME
mylist <- c("JDSP027_stdb", "JDSP089_stdb", "JDSP090_stdb", "JDSP130_stdb", "JDSP138_stdb", "JDSP139_stdb", "JDSP140_stdb", "JDSP141_stdb", "JDSP142_stdb", "JDSP143_stdb", "JDSP144_stdb", "JDSP145_stdb", "JDSP146_stdb", "JDSP147_stdb", "JDSP148_stdb", "JDSP149_stdb", "JDSP150_stdb", "JDSP151_stdb", "JDSP902_stdb", "JDSP909_stdb")

sample_data(FSJphyseq)$label <- NA
sample_data(FSJphyseq)$label[which(sample_data(FSJphyseq)$sample_name %in% mylist)] <- "JDSP"

    

theme(
    plot.title = element_text(color="black", size=20, hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.text.x = element_text(color="black", size=14),
    axis.text.y = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    legend.text = element_text(size=12), legend.title = element_blank())
  
print(p1)

physeq.noncontam.ord <- ordinate(physeq.noncontam, "bray")
p3 = plot_ordination(physeq.noncontam, physeq.noncontam.ord, type="samples", color="sample_type", title="All samples PCA")
print(p1)

# Keep only the most abundant five phyla.
phylum.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:10]
physeq1 = prune_taxa((tax_table(physeq)[, "Phylum"] %in% top5phyla), physeq)

#Plot Ordination
p1 = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord.noncontam, type="taxa", color="Phylum", title="taxa")
print(p1)

physeq.ord.noncontam <- ordinate(physeq2.noncontam1.cleaned, "bray")
plot_ordination(physeq2.noncontam1.cleaned, physeq.ord.noncontam, type="samples", color="sample_type", title="taxa")
print(p1)


physeq1.ord <- ordinate(physeq1, "MDS", "bray", "PCoA")
p1 = plot_ordination(physeq1, physeq1.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

physeq2.noncontam1.cleaned<-prune_samples(FSJsamples,physeq2.noncontam1.cleaned)
physeq.ord <- ordinate(physeq, "MDS", "bray", "PCoA")

physeq.ord.fsj <- ordinate(physeq2.noncontam1.cleaned, "MDS", "bray", "PCoA")
p2 = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord.fsj, type="taxa", color="Phylum") 
print(p2)

physeq1.ord <- ordinate(physeq1, "MDS", "bray", "PCoA")
p1 = plot_ordination(physeq1, physeq1.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

p1 + facet_wrap(~Phylum, 3)

p <- plot_tree(physeq1, color="sample_type", shape="Family", label.tips="Genus", size="abundance")

p2 = plot_ordination(physeq1, physeq1.ord, type="samples", color="physeq1@sam_data[["sample_type"]]") 
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")


##Four Ordination plots
#1.Just OTUs (All Samples)
physeq.ord <- ordinate(physeq2.noncontam1.cleaned, "bray")
p1 = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord, type="taxa", color="Phylum", title="taxa")
print(p1)

p1 + facet_wrap(~Phylum, 3) + scale_fill_manual(values = mycolors)

#2.Samples By Human
p2 = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord, type="samples", color="sample_type", shape="human") 

p2 + geom_polygon(aes(fill=sample_type)) + geom_point(size=5) + ggtitle("samples") 

#2a.Samples All
p2a = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord, type="samples", color="sample_type") 
print(p2a)

p2a + geom_polygon(aes(fill=sample_type)) + geom_point(size=5) + ggtitle("samples") 


#3.Biplot Graphic
p3 = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord, type="biplot", color="sample_type", shape="Phylum", title="biplot")
# Some stuff to modify the automatic shape scale
physeq2.noncontam1.cleaned.shape.names = get_taxa_unique(physeq2.noncontam1.cleaned, "Phylum")
physeq2.noncontam1.cleaned.shape <- 15:(15 + length(physeq2.noncontam1.cleaned.shape.names) - 1)
names(physeq2.noncontam1.cleaned.shape) <- physeq2.noncontam1.cleaned.shape.names
physeq2.noncontam1.cleaned.shape["samples"] <- 16
p3 + scale_shape_manual(values=physeq2.noncontam1.cleaned.shape)

#4.Split Graphic
p4 = plot_ordination(physeq2.noncontam1.cleaned, physeq.ord, type="split", color="Phylum", shape="human", label="sample_type", title="split") 
p4


##Some additional physeq commands
phylum.sum = tapply(taxa_sums(physeqFSJ), tax_table(physeqFSJ)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
physeq1 = prune_taxa((tax_table(physeqFSJ)[, "Phylum"] %in% top5phyla), physeqFSJ)
p<-plot_bar(physeq1, fill = "Phylum")+theme_classic()
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack")


#### FILTER ORDINATION PLOTS PER POPULATION ####
samples<-read.delim("HMP_Metadata_Plus_Skin11Mar22.txt")
ABSsamples<-subset(samples$sample_name, samples$sample_type=="ABS Historic" | samples$sample_type=="ABS Contemporary")
physeqABS<-prune_samples(ABSsamples,physeq2.noncontam1.cleaned)
physeqABS.nozero <- prune_species(taxa_sums(physeqABS)>0, physeqABS)
physeqABS.nozero@sam_data[["sample_type"]]<-as.factor(physeqABS.nozero@sam_data[["sample_type"]])
physeqABS.ord <- ordinate(physeqABS.nozero)
p1 = plot_ordination(physeqABS.nozero, physeqABS.ord, type="samples", color="sample_type")

p1<-p1+ scale_color_manual(values=c("#1b6ca8","#7CBBDE")) + geom_point(size=4) + theme_classic()+
  theme(legend.position = c(0.85, 0.95))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, size=20))+
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust=1, hjust=0), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))
ggsave(p1, file="ABS_PCA.png", width = 7, height = 5, device = "png", dpi=600)
  
  

#### PLE
PLEsamples<-subset(samples$sample_name, samples$sample_type=="PLE Historic" | samples$sample_type=="PLE Contemporary")
physeqPLE<-prune_samples(PLEsamples,physeq2.noncontam1.cleaned)
physeqPLE.nozero <- prune_species(taxa_sums(physeqPLE)>0, physeqPLE)
physeqPLE.nozero@sam_data[["sample_type"]]<-as.factor(physeqPLE.nozero@sam_data[["sample_type"]])
physeqPLE.ord <- ordinate(physeqPLE.nozero)
p2 = plot_ordination(physeqPLE.nozero, physeqPLE.ord, type="samples", color="sample_type")

p2<-p2+ scale_color_manual(values=c("#f37121","#f7ae81")) + geom_point(size=4) + theme_classic()+
  theme(legend.position = c(0.85, 0.95))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, size=20))+
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust=1, hjust=0), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))
p2
ggsave(p2, file="PLE_PCA.png", width = 7, height = 5, device = "png", dpi=600)



#### JDSP
JDSPsamples<-subset(samples$sample_name, samples$sample_type=="JDSP Historic" | samples$sample_type=="JDSP Contemporary")
physeqJDSP<-prune_samples(JDSPsamples,physeq2.noncontam1.cleaned)
physeqJDSP.nozero <- prune_species(taxa_sums(physeqJDSP)>0, physeqJDSP)
physeqJDSP.nozero@sam_data[["sample_type"]]<-as.factor(physeqJDSP.nozero@sam_data[["sample_type"]])
physeqJDSP.ord <- ordinate(physeqJDSP.nozero)
p3 = plot_ordination(physeqJDSP.nozero, physeqJDSP.ord, type="samples", color="sample_type")
p3<-p3+ scale_color_manual(values=c("#c28e04","#f6cd61")) + geom_point(size=4) + theme_classic()+
  theme(legend.position = c(0.85, 0.95))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, size=20))+
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust=1, hjust=0), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))
p3
ggsave(p3, file="JDSP_PCA.png", width = 7, height = 5, device = "png", dpi=600)


#### ONF
ONFsamples<-subset(samples$sample_name, samples$sample_type=="ONF Historic" | samples$sample_type=="ONF Contemporary")
physeqONF<-prune_samples(ONFsamples,physeq2.noncontam1.cleaned)
physeqONF.nozero <- prune_species(taxa_sums(physeqONF)>0, physeqONF)
physeqONF.nozero@sam_data[["sample_type"]]<-as.factor(physeqONF.nozero@sam_data[["sample_type"]])
physeqONF.ord <- ordinate(physeqONF.nozero)
p4 = plot_ordination(physeqONF.nozero, physeqONF.ord, type="samples", color="sample_type")
p4<-p4+ scale_color_manual(values=c("#0c9c5f","#61c3aa")) + geom_point(size=4) + theme_classic()+
  theme(legend.position = c(0.85, 0.95))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, size=20))+
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust=1, hjust=0), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))
p4
ggsave(p4, file="ONF_PCA.png", width = 7, height = 5, device = "png", dpi=600)



### All FSJ -- NOT CLEANED
### NOTE NOTE: Color variable was not found in the available data you provided.No color mapped. ERROR -- YOU JUST NEED ANOTHER COLUMN IN YOUR METADATA

FSJsamples<-read.delim("HMP_Metadata_Plus_Skin11Mar22.txt")
FSJsamples<-FSJsamples$sample_name[1:232]
physeqFSJ.notcleaned<-prune_samples(FSJsamples,physeq)

physeqFSJ.notcleaned.nozero <- prune_species(taxa_sums(physeqFSJ.notcleaned)>0, physeqFSJ.notcleaned)
physeqFSJ.notcleaned.nozero@sam_data[["sample_type"]]<-as.factor(physeqFSJ.notcleaned.nozero@sam_data[["sample_type"]])
physeqFSJ.notcleaned.ord <- ordinate(physeqFSJ.notcleaned.nozero)
p5 = plot_ordination(physeqFSJ.notcleaned.nozero, physeqFSJ.notcleaned.ord, type="samples", color="sample_type")

p5<-p5+ scale_color_manual(values=c("#1b6ca8", "#7CBBDE", "#0c9c5f", "#61C3AA", "#f37121", "#e88f58", "#e09719", "#f6cd61", "#CD4262", "#ff7393")) + geom_point(size=2.5, alpha=0.7) + theme_classic()+
  theme(legend.position = c(0.8, 0.8))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, size=20))+
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust=1, hjust=0), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))
p5

ggsave(p5, file="FSJ_notCleaned_PCA.png", width = 8, height = 7.5, device = "png", dpi=600)




### CLEANED FSJ
physeqFSJ.nozero <- prune_species(taxa_sums(physeqFSJ)>0, physeqFSJ)
physeqFSJ.nozero@sam_data[["sample_type"]]<-as.factor(physeqFSJ.nozero@sam_data[["sample_type"]])
physeqFSJ.nozero.ord <- ordinate(physeqFSJ.nozero)
p6 = plot_ordination(physeqFSJ.nozero, physeqFSJ.nozero.ord, type="samples", color="sample_type")

p6<-p6+ scale_color_manual(values=c("#1b6ca8", "#7CBBDE", "#0c9c5f", "#61C3AA", "#f37121", "#e88f58", "#e09719", "#f6cd61", "#CD4262", "#ff7393")) + geom_point(size=2.5, alpha=0.7) + theme_classic()+
  theme(legend.position = c(0.8, 0.8))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, size=20))+
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14, vjust=1, hjust=0), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=4))+
  theme(plot.margin = margin(15, 22, 15, 22))
p6

ggsave(p6, file="FSJ_Cleaned_PCA.png", width = 8, height = 7.5, device = "png", dpi=600)


##Alpha and Beta Diversity
#plot_richness command: The plot generated by this function will include every sample in physeq, but they can be further grouped on the horizontal axis through the argument to x, and shaded according to the argument to color (see below). You must use untrimmed, non-normalized count data for meaningful results, as many of these estimates are highly dependent on the number of singletons. You can always trim the data later on if needed, just not before using this function
#measure:(Optional). Default is NULL, meaning that all available alpha-diversity measures will be included in plot panels. Alternatively, you can specify one or more measures as a character vector of measure names. Values must be among those supported: c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher").

#Prune OTUs that are not present in any of the samples
physeq2.noncontam1.cleaned.pruned <- prune_taxa(taxa_sums(physeq2.noncontam1.cleaned) > 0, physeq2.noncontam1.cleaned)

#Plot Richness
plot_richness(physeq2.noncontam1.cleaned.pruned)
plot_richness(physeq2.noncontam1.cleaned.pruned, x="SampleType", measures=c("Chao1", "Shannon"))
plot_richness(physeq2.noncontam1.cleaned.pruned, measures=c("Chao1", "Shannon"))


plot_richness(physeq2.noncontam1.cleaned.pruned, measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
plot_richness(physeq2.noncontam1.cleaned.pruned, measures=c("Observed"))

##Organizing Alpha Diversity to a more readable version
#Alpha All Samples
plot_richness(physeq2.noncontam1.cleaned, x="sample_type", measures=c( "Shannon"))

plot_richness(physeq2.noncontam1.cleaned, x="sample_type", measures=c( "Simpson"))

plot_richness(physeq2.noncontam1.cleaned, x="sample_type", measures=c( "Chao1"))

#Alpha By Sample Human
plot_richness(physeq2.noncontam1.cleaned, x="human", color="sample_type", measures=c("Shannon"))

plot_richness(physeq2.noncontam1.cleaned, x="human", color="sample_type", measures=c("Simpson"))

plot_richness(physeq2.noncontam1.cleaned, x="human", color="sample_type", measures=c("Chao1"))


#Alpha By Sample FSJ
plot_richness(physeq2.noncontam1.cleaned, x="FSJ", color="sample_type", measures=c("Shannon"))

plot_richness(physeq2.noncontam1.cleaned, x="FSJ", color="sample_type", measures=c("Simpson"))

plot_richness(physeq2.noncontam1.cleaned, x="FSJ", color="sample_type", measures=c("Chao1"))


#Creating Human and FSJ Sample Types
sampleData(physeq2.noncontam1.cleaned)$human <- getVariable(physeq2.noncontam1.cleaned, "sample_type") %in% c("HMP_oral", "HMP_nasal", "HMP_skin")

sampleData(physeq2.noncontam1.cleaned)$FSJ <- getVariable(physeq2.noncontam1.cleaned, "sample_type") %in% c("ABS Contemporary", "ABS Historic", "JDSP Contemporary", "JDSP Historic", "ONF Contemporary", "ONF Historic", "PLE Contemporary", "PLE Historic", "Savannas")


samp.ord <- ordinate(FSJ.P1,"MDS", "bray", "PCoA")
plot_ordination(FSJ.P1, samp.ord, type="samples", color="sample_type", title="PCoA of FSJ samples after strict contaminant filtering")+ geom_point(size=3, alpha=0.9)


#Types of Alpha Diversity measures
#Chao1:–	estimate	diversity	from	abundance	data	(importance	of	rare	OTUs)
#ACE: involves an arbitrary	abundance	threshold to label Sabun	as the number of abundant taxa, Srare	as the	number of rare	taxa; The expression basically inflates the number of rare taxa	and inflates again	the	number of taxa with	abundance	1.
#Shannon: How difficult it is to predict the identity of a randomly chosen individual.
#Simpson: The probability that two randomly chosen individuals are the same species.
#Inverse Simpson: Assuming a theoretically community where all species were equally abundant, this would be the number of species needed to have the same Simpson index value for the community being analyzed.
