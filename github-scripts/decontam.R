## REMOVING HMP TAXA CONTAMINANTS FROM OTU TABLE - Phyloseq and Decontam
## QC sequencing reads and taxa
## Dec 2022

setwd("~/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/")
#devtools::install_github("vmikk/metagMisc")
#install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

#devtools::install_github("david-barnett/microViz")

library(qiime2R)
library(decontam)
library(biomformat);packageVersion("biomformat")
library(ggplot2)
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library(phyloseq)
library(microbiome)
library(knitr)
library(remotes)
library(metagMisc)
library("data.table"); packageVersion("data.table")
library(microViz)

source("~/Box Sync/PhD/Projects/Kraken_Pathogen/removing_HMP/phyloseq_funcs.R", local = TRUE) ### load in functions for analysis

## Outputs are from Qiime2
physeq<-qza_to_phyloseq(
  features= "HMP-table-no-hsapiens-exact.qza",
  taxonomy = "HMP-feature-tax-table-2022.qza",
  metadata = "FSJ_HMP_Metadata_table_2022.txt") #5736 unfiltered taxa


##### Filtering taxa by minimum abundance  or frequency across samples!
# Note: in 16s data: We usually take a quick first look at the library sizes (i.e. the number of reads) in each sample. 
# Note: our data is metagenomic!

# OTUs that are found in at least 5% of samples
# phyloseq_filter_prevalence(physeq, prev.trh = 0.05, abund.trh = NULL)

# prev.trh = Prevalence threshold (default, 0.05 = 5% of samples)
# threshold_condition = Indicates type of prevalence and abundance conditions, can be "OR" (default) or "AND"
# abund.type = Character string indicating which type of OTU abundance to take into account for filtering ("total", "mean", or "median")
physeq <- phyloseq_filter_prevalence(physeq, prev.trh = 0, abund.trh = 5, threshold_condition = "AND", abund.type = "total")  # 4447 taxa left
# total OTU abundance is >= 5 reads it'll be preserved too

# Include only taxa with more than 10 reads (on average) in at least 10% samples
# phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.1, abund.trh = 10, abund.type = "mean", threshold_condition = "AND")  # 4447 taxa left


#########################  EXPLORE TOTAL KRAKEN READ DEPTH FOUND ############################
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

depth.plot <- pSeqDepth + facet_wrap(~Sample_or_Control) + ggtitle("Total Read Counts")
depth.plot


#### Taxa total counts histogram
# 1. How many singleton taxa are there (OTUs that occur in just one sample, one time)?
# 2. How many doubletons are there (OTU that occurs just twice)?
# 3. Create a histogram of the total counts of each OTU.
# 4. Calculate the cumulative sum of OTUs that would be filtered at every possible value of such a threshold, from zero to the most-observed OTU. This one is tricky. Feel free to glance at the answers. I used some data.table magic to make this easier
# 5. Plot the cumulative sum of OTUs against the total counts using ggplot2 commands to make a scatter plot, and save this object as pCumSum, then “print” it to the terminal to render a graphic. What behavior do you see in the data? Are most of the OTUs in the table rare? Where would you set the threshold?
# 6. To help clarify, zoom-in on the region between zero and 100 total counts, by “adding” the following layer to pCumSum: pCumSum + xlim(0, 100) Now where would you set the threshold?
  
tdt = data.table(tax_table(physeq),
                 TotalCounts = taxa_sums(physeq),
                 OTU = taxa_names(physeq))

taxaCount<-ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Taxa Counts")

taxaCount

# How many singletons?
tdt[(TotalCounts <= 0), .N]

# microbes with observations < 5 (should be zero because we set a minimum # of reads)
tdt[(TotalCounts < 5), .N]

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum
pCumSum + xlim(0, 100)


# Taxa prevalence histogram
mdt = fast_melt(physeq)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID] #function “melts” the OTU table into a long format with three main columns: SampleID, TaxaID, and count

ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Histogram of Taxa Prevalence")

# singletons?
prevdt[(Prevalence <= 0), .N]
prevdt[(Prevalence <= 2), .N]
prevdt[(Prevalence <= 5), .N] #1463 will be filtered out

# Prevalence vs. Total Count Scatter plot
ggplot(prevdt, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 2, alpha = 0.75) + 
  scale_y_log10()


# Summarize info from your data
summarize_phyloseq(physeq)
meta <- meta(physeq)
taxonomy <- tax_table(physeq)

# Absolute abundances
otu.absolute <- abundances(physeq)

# Relative abundances
otu.relative <- abundances(physeq, "compositional")

# Total read counts:
reads_sample <- readcount(physeq)
# check for first 5 samples
reads_sample[1:5]

# We will now plot the total number of counts per sample.
sample_sums(physeq) #get read count per sample
sort(sample_sums(physeq))
hist(sample_sums(physeq), main="Histogram: Read Counts", xlab="Total Reads", las=1, breaks=12)
summary(sample_sums(physeq))

metadata <- data.frame(sample_data(physeq)) # add this total read count to your metadata
head(metadata)

metadata$total_reads <- sample_sums(physeq)



################################################### REMOVE HMP SPECIES ##########################################################
# Decontam Removal -- Method: Frequency
# The first contaminant identification method we’ll use is the “frequency” method. In this method, the distribution of the frequency of each sequence feature as a function of the input DNA concentration is used to identify contaminants.
# contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
# head(contamdf.freq)

# Decontam Removal -- Method:Prevalence
# The second contaminant identification method we’ll use is the “prevalence” method. In this method, the prevalence (presence/absence across samples) of each sequence feature in true positive samples is compared to the prevalence in negative controls to identify contaminants.

# In our phyloseq object, "Sample_or_Control" is the sample variable that holds the negative control information. We’ll summarize that data as a logical variable, with TRUE for control samples, as that is the form required by isContaminant.

sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant) ## Default is threshold of 0.1

# FALSE  TRUE 
# 3015   1432
# filter out "true contaminant"

# Make phyloseq object of presence-absence in negative controls and true samples
physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$Sample_or_Control == "Control Sample", physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$Sample_or_Control == "True Sample", physeq.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                    contaminant=contamdf.prev$contaminant)

p <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative HMP Controls)") + ylab("Prevalence (True FSJ Samples)")+ ggtitle("Contaminant Threshold 0.1") +
  scale_color_manual(values=c("#5F9DF7", "#1746A2"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme( # remove the vertical grid lines 
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.1, color="black" )) +
  theme(legend.title = element_text(hjust=0.5, size=14))+
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x=element_text(size = 14), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=5))+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(plot.title = element_text(color="black", size=16),)+
  theme(plot.margin = margin(12, 12, 10, 22))
p
ggsave(p, filename = "contam_threshold0.1_scatter.png", width=8.5, height=6, units = "in", device="png", dpi=800)


#threshold 0.5 plot
contamdf.prev05 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

# FALSE  TRUE 
# 2510  1937

df.pa2 <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                     contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa2, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+ ggtitle("Contaminant Threshold 0.5")


#threshold 1 plot -- MOST STRINGENT
contamdf.prev1 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=1)
table(contamdf.prev1$contaminant)

#FALSE  TRUE 
#1527  2920

df.pa3 <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                     contaminant=contamdf.prev1$contaminant)
ggplot(data=df.pa3, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+ ggtitle("Contaminant Threshold 1.0")



# Make phyloseq object of presence-absence in negative controls and true samples
physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$Sample_or_Control == "Control Sample", physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$Sample_or_Control == "True Sample", physeq.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg),
                    contaminant=contamdf.prev1$contaminant)

p <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative HMP Controls)") + ylab("Prevalence (True FSJ Samples)")+ ggtitle("Contaminant Threshold 1") +
  scale_color_manual(values=c("#5F9DF7", "#1746A2"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme( # remove the vertical grid lines 
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line( size=.1, color="black" )) +
  theme(legend.title = element_text(hjust=0.5, size=14))+
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text.x=element_text(size = 14), axis.text.y=element_text(size = 14))+
  theme(axis.title.y = element_text(vjust=5))+
  theme(axis.title.x = element_text(vjust=-0.5))+
  theme(plot.title = element_text(color="black", size=16),)+
  theme(plot.margin = margin(12, 12, 10, 22))
p
ggsave(p, filename = "contam_threshold1_scatter.png", width=8.5, height=6, units = "in", device="png", dpi=800)



########################### PRUNE -- Command to filter out contaminants
physeq.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
physeq.noncontam

physeq.noncontam.P1 <- prune_taxa(!contamdf.prev1$contaminant, physeq)
physeq.noncontam.P1

## Look at actual Contaminants
Contam.list <- subset(contamdf.prev,contamdf.prev$contaminant =="TRUE")
Contam.list<-tibble::rownames_to_column(Contam.list, "TaxonID")
Contam.list

## Determine what species each taxonID corresponds to in Contaminants
library(tibble)
physeq_taxa_df<-as.data.frame(physeq@tax_table@.Data)
df <- tibble::rownames_to_column(physeq_taxa_df, "TaxonID")
contam_taxon<-subset(df,df$TaxonID %in%Contam.list$TaxonID)

contam_taxon_edit<-data.frame(table(contam_taxon$Kingdom))

## Contamination Visualizations:
## by Kingdom
bp<-ggplot(contam_taxon_edit, aes(x = "", y = Freq, fill = Var1))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
pie

pie + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=18, hjust = 0.5, vjust = 0.5),
    legend.text=element_text(size=14),
    legend.title=element_blank()
  )

pie + scale_fill_grey() +  blank_theme +
  ggtitle("Proportion of Taxa Removed\n(Kingdom)")+
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


## By Genus
contam_taxon_Genus<-data.frame(table(contam_taxon$Genus))
contam_taxon_Genus <- contam_taxon_Genus[order(-contam_taxon_Genus$Freq),]
contam_taxon_Genus_top <- contam_taxon_Genus[1:50,]

contam_taxon_Genus_top$Var1 <- factor(contam_taxon_Genus_top$Var1, levels = contam_taxon_Genus$Var1[1:50])

genusout<-ggplot(contam_taxon_Genus_top, aes(x = Var1, y = Freq))+
  geom_bar(stat = "identity", fill='steelblue') +
  theme(axis.text.x = element_text(angle = 90))+ 
  ggtitle("Top 50 contaminants that were filtered out\n(Genus)") +
  xlab("Taxonomic name (Genus)") + ylab("Count")+ theme(
    plot.title = element_text(color="black", size=20, hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.text.x = element_text(color="black", size=14),
    axis.text.y = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14))
### May want to color this by Kingdom later on...

ggsave(genusout, filename = "../removing_HMP/50topGenus_filtered.pdf", width = 16, height = 9, units = "in", device="pdf", dpi=800)


## By Species
contam_taxon_spps<-data.frame(table(contam_taxon$Species))
contam_taxon_spps <- contam_taxon_spps[order(-contam_taxon_spps$Freq),]
contam_taxon_spps_top <- contam_taxon_spps[1:50,]

contam_taxon_spps_top$Var1 <- factor(contam_taxon_spps_top$Var1, levels = contam_taxon_spps$Var1[1:50])

ggplot(contam_taxon_spps_top, aes(x = Var1, y = Freq))+
  geom_bar(stat = "identity", fill='steelblue') +
  theme(axis.text.x = element_text(angle = 90))+ 
  ggtitle("Top 50 contaminants that were filtered out\n(Species)") +
  xlab("Taxonomic name") + ylab("Count")+ theme(
    plot.title = element_text(color="black", size=20, hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.text.x = element_text(color="black", size=14),
    axis.text.y = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14))
### May want to color this by Kingdom later on...


contam_taxon$Genus[contam_taxon$Species %in% contam_taxon_spps_top$Var1]


######################################## remove HMP samples before exporting .biom table #########################################
FSJphyseq_cleaned = ps_filter(physeq.noncontam, Sample_or_Control != "Control Sample", .keep_all_taxa = F)# filter out the HMP data
### NOTE: human reads have already been filtered out in the original Qiime2 object 

## check that there are no taxa without any reads
rowSums(FSJphyseq_cleaned@otu_table)


FSJ.P1 = ps_filter(physeq.noncontam.P1, Sample_or_Control != "Control Sample", .keep_all_taxa = F)
rowSums(FSJ.P1@otu_table)


################################### Save as R object #####################################
saveRDS(FSJphyseq_cleaned, file = "noHMP.noSapien.rds")
saveRDS(physeq, file="raw_HMP_FSJ_physeq.rds")

################################### How to extract phyloseq OTU #######################################
# Prepare and Export Taxonomy, OTU Table, and Metadata
# https://forum.qiime2.org/t/importing-dada2-and-phyloseq-objects-to-qiime-2/4683

# Export taxonomy table as "tax.txt"

tax<-as(tax_table(FSJphyseq_cleaned),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "tax_decontam2022.txt", quote=FALSE, col.names=FALSE, sep="\t")

# Export feature/OTU table
# As a biom file

library(biomformat);packageVersion("biomformat")
## [1] ‘1.6.0’

otu<-t(as(otu_table(FSJphyseq_cleaned),"matrix")) # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
#otu<-as(otu_table(GlobalPatterns),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_decontam2022.biom") ## Make into .biom table for Qiime2 and for SourceTracker!

# As a text file
#write.table(t(seqtab), "seqtab.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
#or from the phyloseq object, 't' to transform if taxa_are_rows=FALSE, no 't' if taxa_are_rows=TRUE
#write.table(t(otu_table(ps), "seqtab.txt",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Export metadata (if you have a properly formatted metadata file that you imported in your phyloseq pipeline, you can skip this step and just use that text file directly in QIIME 2)

write.table(sample_data(newPhyloObject),"sample-metadata-decontam2022.txt", sep="\t", row.names=T, col.names=TRUE, quote=FALSE)
# NOTE: have to manipulate in excel to add "sample_name" as the first column header