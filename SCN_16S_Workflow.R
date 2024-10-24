# Emily A. Green, PhD. 
# Code for 16S data, includes metaphlan taxonomy for cyst metagenomes and kraken2 taxonomy for cyst and soil metagenomes 

#16S Cyst samples that are not present in the cyst metagenome were removed
#All S1 samples were removed

#This code reflects being run in the /workdir/eag252/LTR_16S_Metagenome/ directory.
#If being run in the "storage" directory, file paths will likely need to be updated. 

## dada2 processing-----
library(dada2)
path <- "/workdir/eag252/LTR_16S_Metagenome/"
list.files(path)

#### Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[40:60])
plotQualityProfile(fnRs[40:60])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#default parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "/workdir/eag252/LTR_16S_Metagenome/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

#Exports file that shows read counts after each of the dada2 processing steps
write.csv(track, "all16scounts.csv")

#Output the taxa, taxa counts, and taxa asv sequences
write.csv(taxa.print, "taxaLTR16S.csv")
write.csv(seqtab.nochim, "ASVLTR16S.csv")

#Load packages----
library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)
library(viridisLite)
library(viridis)
library(ape)
library(cowplot)
theme_set(theme_bw())
library(NetCoMi)

#Read in metadata----
meta <- read.csv("/workdir/eag252/LTR_16S_Metagenome/LTR16s_S3ARSMetadata.csv",header = TRUE, row.names = 1)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(meta), tax_table(taxa))
sample_data(ps)
ps


#Filter taxa----
#Confirm all reads are identified to a Kingdom
sort(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))

#Remove reads that are not identified to "Bacteria"
ps1 <- subset_taxa(ps,Kingdom == "Bacteria")
get_taxa_unique(ps1, taxonomic.rank = "Kingdom")
ps1

sort(get_taxa_unique(ps1, taxonomic.rank = "Phylum"))
sort(get_taxa_unique(ps1, taxonomic.rank = "Family"))
sort(get_taxa_unique(ps1, taxonomic.rank = "Order"))

#Remove mitochondria & chloroplast reads
ps2 <- subset_taxa(ps1, !Family == "Mitochondria")
sort(get_taxa_unique(ps2, taxonomic.rank = "Family"))
ps2 <- subset_taxa(ps2, !Order == "Chloroplast")
ps2

#Decontam----
#Checks for potential contaminates that are in the negative controls & present in the samples
library(decontam)
sample_variables(ps2)
df <- as.data.frame(sample_data(ps2))
df$LibrarySize <- sample_sums(ps2)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
sample_data(ps2)$is.neg <- sample_data(ps2)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps2, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant)

sample_sums(ps2)
View(sample_sums(ps2))

#Remove samples with less than 100 reads----
ps3 <- prune_samples(sample_sums(ps2) > 100, ps2)
View(sample_sums(ps3))


#Check Beta Diversity with Unifrac. Low read samples removed ----
ps3t <- ps3

#Unifrac
ps3t.tree <- rtree(ntaxa(ps3t), rooted = TRUE, tip.label = taxa_names(ps3t))
#plot(ps3t.tree)
ps3t <- merge_phyloseq(ps3t, ps3t.tree)

#Unweighted
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=F)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "Sample_or_Control", shape = "Sample_Type") + theme(aspect.ratio=1)

#statistics
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~Sample_Type*Sample_or_Control, ps3t.tree.data)

#Weighted
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=T)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "Sample_or_Control", shape = "Sample_Type") + theme(aspect.ratio=1)

#statistics
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~Sample_Type*Sample_or_Control, ps3t.tree.data)

#Make new phyloseq object. Remove all controls & mock communities---- 
#ps4 has no negative controls or mock communities
#Remove controls
ps4 <- subset_samples(ps3, Sample_or_Control == "Sample")
View(sample_data(ps4))

#Beta diversity with Unifrac- All samples with no negative controls or mock communities----
ps4t <-ps4 

#Unifrac
ps4t.tree <- rtree(ntaxa(ps4t), rooted = TRUE, tip.label = taxa_names(ps4t))
#plot(ps4t.tree)
ps4t <- merge_phyloseq(ps4t, ps4t.tree)

#Unweighted
ps4tunifrac_dist <- phyloseq::distance(ps4t, method="unifrac", weighted=F)
ordination = ordinate(ps4t, method="PCoA", distance=ps4tunifrac_dist)
plot_ordination(ps4t, ordination, color = "Sample_Type", shape = "Season") + theme(aspect.ratio=1) + geom_point(size = 2)

#statistics
#permanova
ps4t.tree.data <- data.frame(sample_data(ps4t))
adonis2(ps4tunifrac_dist ~Sample_Type*Season*Plot, ps4t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(ps4tunifrac_dist, ps4t.tree.data$Sample_Type, p.adjust.m = "holm", perm = 9999)

pairwiseAdonis::pairwise.adonis(ps4tunifrac_dist, ps4t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
ps4tu.bdisper <- betadisper(ps4tunifrac_dist, ps4t.tree.data$Sample_Type, "centroid", bias.adjust = TRUE)
anova(ps4tu.bdisper)
permutest(ps4tu.bdisper, pairwise = TRUE)


#Weighted
ps4tunifrac_dist <- phyloseq::distance(ps4t, method="unifrac", weighted=T)
ordination = ordinate(ps4t, method="PCoA", distance=ps4tunifrac_dist)
plot_ordination(ps4t, ordination, color = "Sample_Type", shape = "Season") + theme(aspect.ratio=1) + geom_point(size = 2)

#statistics
#permanova
ps4t.tree.data <- data.frame(sample_data(ps4t))
adonis2(ps4tunifrac_dist ~Sample_Type*Season*Plot, ps4t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(ps4tunifrac_dist, ps4t.tree.data$Sample_Type, p.adjust.m = "holm", perm = 9999)

pairwiseAdonis::pairwise.adonis(ps4tunifrac_dist, ps4t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
ps4tu.bdisper <- betadisper(ps4tunifrac_dist, ps4t.tree.data$Sample_Type, "centroid", bias.adjust = TRUE)
anova(ps4tu.bdisper)
permutest(ps4tu.bdisper, pairwise = TRUE)


# Split samples into separate phyloseq objects----
pssoil <- subset_samples(ps4, SampleType == "Soil")

psbulk <- subset_samples(ps4, SoilType == "Bulk soil")

psrhiz <- subset_samples(ps4, SoilType == "Rhizosphere")

pscyst <- subset_samples(ps4, SampleType == "Cyst")

psroot <- subset_samples(ps4, SampleType == "Roots")


#Remove taxa with 0 reads from each of the individual sample type phyloseq objects----
#6782 taxa
pssoil <- prune_taxa(taxa_sums(pssoil) > 0, pssoil)

psbulk <- prune_taxa(taxa_sums(psbulk) > 0, psbulk)
View(tax_table(psbulk))

psrhiz <- prune_taxa(taxa_sums(psrhiz) > 0, psrhiz)
View(tax_table(psrhiz))

pscyst <- prune_taxa(taxa_sums(pscyst) > 0, pscyst)
View(tax_table(pscyst))

psroot <- prune_taxa(taxa_sums(psroot) > 0, psroot)
View(tax_table(psroot))

#Alpha Diversity. Observed, Shannon & Simpson----
#Need to split the shannon and simpson plots to view the box plots without the individual dots. 
#https://stackoverflow.com/questions/73544659/geom-points-are-not-placed-on-the-boxplot

#ps4 object- all samples
ps4alpha <- estimate_richness(ps4, measures = c("Observed", "Shannon", "Simpson"))
ps4alpha$Sample_Type <- sample_data(ps4)$Sample_Type
ps4alpha$Plot <- sample_data(ps4)$Plot
ps4alpha$Season <- sample_data(ps4)$Season
ps4alpha$SoilType <- sample_data(ps4)$SoilType


ggplot(ps4alpha, aes(x = Plot, y = Observed, color = Season)) + geom_boxplot()+ facet_grid(~Sample_Type, scales = "free", space = "free" )
#statistics 
ps4alpha.sampletype <- aov(Observed ~Sample_Type, ps4alpha)
anova(ps4alpha.sampletype)
ps4alpha.season <- aov(Observed ~Season, ps4alpha)
anova(ps4alpha.season)
ps4alpha.plot <- aov(Observed ~Plot, ps4alpha)
anova(ps4alpha.plot)


ggplot(ps4alpha, aes(x = Plot, y = Shannon, color = Season)) + geom_boxplot()+ facet_grid(~Sample_Type, scales = "free", space = "free" )
#statistics 
ps4alpha.sampletype <- aov(Shannon ~Sample_Type, ps4alpha)
anova(ps4alpha.sampletype)
ps4alpha.season <- aov(Shannon ~Season, ps4alpha)
anova(ps4alpha.season)
ps4alpha.plot <- aov(Shannon ~Plot, ps4alpha)
anova(ps4alpha.plot)


ggplot(ps4alpha, aes(x = Plot, y = Simpson, color = Season)) + geom_boxplot()+ facet_grid(~Sample_Type, scales = "free", space = "free" )
#statistics 
ps4alpha.sampletype <- aov(Simpson ~Sample_Type, ps4alpha)
anova(ps4alpha.sampletype)
ps4alpha.season <- aov(Simpson ~Season, ps4alpha)
anova(ps4alpha.season)
ps4alpha.plot <- aov(Simpson ~Plot, ps4alpha)
anova(ps4alpha.plot)


#Do a dunn's test when the anova is significant (aka the means are not equal)
library(dunn.test)
dunn.test(ps4alpha$Observed, g = ps4alpha$Sample_Type, method = "bonferroni")

dunn.test(ps4alpha$Shannon, g = ps4alpha$Sample_Type, method = "bonferroni")

dunn.test(ps4alpha$Simpson, g = ps4alpha$Sample_Type, method = "bonferroni")



#Soil
#plot_richness(pssoil, x="Plot", measures = c("Shannon", "Simpson"), color = "SoilType", shape = "Season") 
soilalpha <- estimate_richness(pssoil, measures = c("Observed", "Shannon", "Simpson"))
soilalpha$Plot <- sample_data(pssoil)$Plot
soilalpha$SoilType <- sample_data(pssoil)$SoilType
soilalpha$Season <- sample_data(pssoil)$Season

ggplot(soilalpha, aes(x = Plot, y = Observed, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
#statistics 
soilalpha.soiltype <- aov(Observed ~SoilType, soilalpha)
anova(soilalpha.soiltype)
soilalpha.season <- aov(Observed ~Season, soilalpha)
anova(soilalpha.season)
soilalpha.plot <- aov(Observed ~Plot, soilalpha)
anova(soilalpha.plot)

#dunns test
dunn.test(soilalpha$Observed, g = soilalpha$Sample_Type, method = "bonferroni")

dunn.test(soilalpha$Observed, g = soilalpha$Season, method = "bonferroni")


ggplot(soilalpha, aes(x = Plot, y = Shannon, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
#statistics 
soilalpha.soiltype <- aov(Shannon ~SoilType, soilalpha)
anova(soilalpha.soiltype)
soilalpha.season <- aov(Shannon ~Season, soilalpha)
anova(soilalpha.season)
soilalpha.plot <- aov(Shannon ~Plot, soilalpha)
anova(soilalpha.plot)


ggplot(soilalpha, aes(x = Plot, y = Simpson, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
#statistics 
soilalpha.soiltype <- aov(Simpson ~SoilType, soilalpha)
anova(soilalpha.soiltype)
soilalpha.season <- aov(Simpson ~Season, soilalpha)
anova(soilalpha.season)
soilalpha.plot <- aov(Simpson ~Plot, soilalpha)
anova(soilalpha.plot)


#Cyst
#plot_richness(pscyst, x="Plot", measures = c("Shannon", "Simpson"), color = "Season")
cystalpha <- estimate_richness(pscyst,measures = c("Observed", "Shannon", "Simpson"))
cystalpha$Plot <- sample_data(pscyst)$Plot
cystalpha$Season <- sample_data(pscyst)$Season

ggplot(cystalpha, aes(x= Plot, y = Observed, color = Season)) + geom_boxplot()
#statistics 
cystalpha.plot <- aov(Observed ~Plot, cystalpha)
anova(cystalpha.plot)
cystalpha.season <- aov(Observed ~Season, cystalpha)
anova(cystalpha.season)

#dunns test
dunn.test(cystalpha$Observed, g = cystalpha$Plot, method = "bonferroni")


ggplot(cystalpha, aes(x= Plot, y = Shannon, color = Season)) + geom_boxplot()
#statistics 
cystalpha.plot <- aov(Shannon ~Plot, cystalpha)
anova(cystalpha.plot)
cystalpha.season <- aov(Shannon ~Season, cystalpha)
anova(cystalpha.season)

ggplot(cystalpha, aes(x= Plot, y = Simpson, color = Season)) + geom_boxplot()
#statistics 
cystalpha.plot <- aov(Simpson ~Plot, cystalpha)
anova(cystalpha.plot)
cystalpha.season <- aov(Simpson ~Season, cystalpha)
anova(cystalpha.season)


#Root
#plot_richness(psroot, x="Plot", measures = c("Shannon", "Simpson"), color = "Season")
rootalpha <- estimate_richness(psroot, measures = c("Observed", "Shannon", "Simpson"))
rootalpha$Plot <- sample_data(psroot)$Plot

ggplot(rootalpha, aes(x= Plot, y = Observed)) + geom_boxplot()
#statistics 
rootalpha.plot <- aov(Observed ~Plot, rootalpha)
anova(rootalpha.plot)

ggplot(rootalpha, aes(x= Plot, y = Shannon)) + geom_boxplot()
rootalpha.plot <- aov(Shannon ~Plot, rootalpha)
anova(rootalpha.plot)

ggplot(rootalpha, aes(x= Plot, y = Simpson)) + geom_boxplot()
#statistics 
rootalpha.plot <- aov(Simpson ~Plot, rootalpha)
anova(rootalpha.plot)

#Beta Diversity using Unifrac. ----

#Create new phyloseq object for adding phylogenetic tree
#Certain downstream analyses won't work with a phylogenetic tree in the phyloseq object.
#I just create a new phyloseq object for this reason. 
pssoil.t <- pssoil
pscyst.t <- pscyst
psroot.t <- psroot
psrhiz.t <- psrhiz
psbulk.t <- psbulk

psbulkmid <- subset_samples(psbulk, Season == "Mid")
psbulkmid.t <- psbulkmid
psbulkfall <- subset_samples(psbulk, Season == "Fall")
psbulkfall.t <- psbulkfall

pscystmid <- subset_samples(pscyst, Season == "Mid")
pscystmid.t <- pscystmid
pscystfall <- subset_samples(pscyst, Season == "Fall")
pscystfall.t <- pscystfall


#Beta diversity- Unifrac Soil (bulk soil & rhizosphere) ----
pssoil.t.tree <- rtree(ntaxa(pssoil.t), rooted = TRUE, tip.label = taxa_names(pssoil.t))
#plot(pssoil.t.tree)
pssoil.t <- merge_phyloseq(pssoil.t, pssoil.t.tree)

#Unweighted Soil
pssoiltunifrac_dist <- phyloseq::distance(pssoil.t, method="unifrac", weighted=F)
ordination = ordinate(pssoil.t, method="PCoA", distance=pssoiltunifrac_dist)
plot_ordination(pssoil.t, ordination, color = "SoilType", shape = "Season", label = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.t, ordination, color = "Plot", shape = "SoilType", label = "Season") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pssoil.t.tree.data <- data.frame(sample_data(pssoil.t))
adonis2(pssoiltunifrac_dist ~SoilType*Plot*Season, pssoil.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pssoiltunifrac_dist, pssoil.t.tree.data$SoilType, p.adjust.m = "holm", perm = 9999)
pairwiseAdonis::pairwise.adonis(pssoiltunifrac_dist, pssoil.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
pssoiltu.bdisper <- betadisper(pssoiltunifrac_dist, pssoil.t.tree.data$SoilType, "centroid", bias.adjust = TRUE)
anova(pssoiltu.bdisper)
permutest(pssoiltu.bdisper, pairwise = TRUE)

#Weighted Soil
pssoiltunifrac_dist <- phyloseq::distance(pssoil.t, method="unifrac", weighted=T)
ordination = ordinate(pssoil.t, method="PCoA", distance=pssoiltunifrac_dist)
plot_ordination(pssoil.t, ordination, color = "SoilType", shape = "Season", label = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.t, ordination, color = "Plot", shape = "SoilType", label = "Season") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pssoil.t.tree.data <- data.frame(sample_data(pssoil.t))
adonis2(pssoiltunifrac_dist ~SoilType*Plot*Season, pssoil.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pssoiltunifrac_dist, pssoil.t.tree.data$SoilType, p.adjust.m = "holm", perm = 9999)
pairwiseAdonis::pairwise.adonis(pssoiltunifrac_dist, pssoil.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

pssoiltu.bdisper <- betadisper(pssoiltunifrac_dist, pssoil.t.tree.data$SoilType, "centroid", bias.adjust = TRUE)
anova(pssoiltu.bdisper)
permutest(pssoiltu.bdisper, pairwise = TRUE)

#Beta diversity- Unifrac mid season soil(bulk soil & rhizosphere) ---- 

pssoil.mid <- subset_samples(pssoil, Season == "Mid")
pssoil.mid.t <- pssoil.mid

pssoil.mid.t.tree <- rtree(ntaxa(pssoil.mid.t), rooted = TRUE, tip.label = taxa_names(pssoil.mid.t))
#plot(pssoil.mid.t.tree)
pssoil.mid.t <- merge_phyloseq(pssoil.mid.t, pssoil.mid.t.tree)

#Unweighted Soil
pssoil.midtunifrac_dist <- phyloseq::distance(pssoil.mid.t, method="unifrac", weighted=F)
ordination = ordinate(pssoil.mid.t, method="PCoA", distance=pssoil.midtunifrac_dist)
plot_ordination(pssoil.mid.t, ordination, color = "SoilType", shape = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.mid.t, ordination, color = "Plot", shape = "SoilType") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pssoil.mid.t.tree.data <- data.frame(sample_data(pssoil.mid.t))
adonis2(pssoil.midtunifrac_dist ~SoilType*Plot, pssoil.mid.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$SoilType, p.adjust.m = "holm", perm = 9999)
pairwiseAdonis::pairwise.adonis(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
pssoiltmidu.bdisper <- betadisper(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$SoilType, "centroid", bias.adjust = TRUE)
anova(pssoiltmidu.bdisper)
permutest(pssoiltmidu.bdisper, pairwise = TRUE)

#Weighted Soil
pssoil.midtunifrac_dist <- phyloseq::distance(pssoil.mid.t, method="unifrac", weighted=T)
ordination = ordinate(pssoil.mid.t, method="PCoA", distance=pssoil.midtunifrac_dist)
plot_ordination(pssoil.mid.t, ordination, color = "SoilType", shape = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.mid.t, ordination, color = "Plot", shape = "SoilType") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pssoil.mid.t.tree.data <- data.frame(sample_data(pssoil.mid.t))
adonis2(pssoil.midtunifrac_dist ~SoilType*Plot, pssoil.mid.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$SoilType, p.adjust.m = "holm", perm = 9999)
pairwiseAdonis::pairwise.adonis(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
pssoilmidtu.bdisper <- betadisper(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$SoilType, "centroid", bias.adjust = TRUE)
anova(pssoiltmidu.bdisper)
permutest(pssoiltmidu.bdisper, pairwise = TRUE)

#Beta diversity- Unifrac Rhizosphere ----
psrhiz.t.tree <- rtree(ntaxa(psrhiz.t), rooted = TRUE, tip.label = taxa_names(psrhiz.t))
#plot(psrhiz.t.tree)
psrhiz.t <- merge_phyloseq(psrhiz.t, psrhiz.t.tree)

#Unweighted Rhizosphere
psrhiztunifrac_dist <- phyloseq::distance(psrhiz.t, method="unifrac", weighted=F)
ordination = ordinate(psrhiz.t, method="PCoA", distance=psrhiztunifrac_dist)
plot_ordination(psrhiz.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psrhiz.t.tree.data <- data.frame(sample_data(psrhiz.t))
adonis2(psrhiztunifrac_dist ~Plot, psrhiz.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psrhiztunifrac_dist, psrhiz.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
psrhiztu.bdisper <- betadisper(psrhiztunifrac_dist, psrhiz.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psrhiztu.bdisper)
permutest(psrhiztu.bdisper, pairwise = TRUE)



#Weighted Rhizosphere
psrhiztunifrac_dist <- phyloseq::distance(psrhiz.t, method="unifrac", weighted=T)
ordination = ordinate(psrhiz.t, method="PCoA", distance=psrhiztunifrac_dist)
plot_ordination(psrhiz.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psrhiz.t.tree.data <- data.frame(sample_data(psrhiz.t))
adonis2(psrhiztunifrac_dist ~Plot, psrhiz.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psrhiztunifrac_dist, psrhiz.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

psrhiztu.bdisper <- betadisper(psrhiztunifrac_dist, psrhiz.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psrhiztu.bdisper)
permutest(psrhiztu.bdisper, pairwise = TRUE)



#Beta diversity- Unifrac Bulk Soil ----
psbulk.t.tree <- rtree(ntaxa(psbulk.t), rooted = TRUE, tip.label = taxa_names(psbulk.t))
#plot(psbulk.t.tree)
psbulk.t <- merge_phyloseq(psbulk.t, psbulk.t.tree)

#Unweighted Bulk
psbulktunifrac_dist <- phyloseq::distance(psbulk.t, method="unifrac", weighted=F)
ordination = ordinate(psbulk.t, method="PCoA", distance=psbulktunifrac_dist)
plot_ordination(psbulk.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psbulk.t.tree.data <- data.frame(sample_data(psbulk.t))
adonis2(psbulktunifrac_dist ~Plot*Season, psbulk.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psbulktunifrac_dist, psbulk.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
psbulktu.bdisper <- betadisper(psbulktunifrac_dist, psbulk.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulktu.bdisper)
permutest(psbulktu.bdisper, pairwise = TRUE)


#Weighted Bulk
psbulktunifrac_dist <- phyloseq::distance(psbulk.t, method="unifrac", weighted=T)
ordination = ordinate(psbulk.t, method="PCoA", distance=psbulktunifrac_dist)
plot_ordination(psbulk.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psbulk.t.tree.data <- data.frame(sample_data(psbulk.t))
adonis2(psbulktunifrac_dist ~Plot*Season, psbulk.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psbulktunifrac_dist, psbulk.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
psbulktu.bdisper <- betadisper(psbulktunifrac_dist, psbulk.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulktu.bdisper)
permutest(psbulktu.bdisper, pairwise = TRUE)



#Beta diversity- Unifrac Bulk Soil Mid season samples ----
psbulkmid.t.tree <- rtree(ntaxa(psbulkmid.t), rooted = TRUE, tip.label = taxa_names(psbulkmid.t))
#plot(psbulkmid.t.tree)
psbulkmid.t <- merge_phyloseq(psbulkmid.t, psbulkmid.t.tree)

#Unweighted Bulk Mid
psbulkmidtunifrac_dist <- phyloseq::distance(psbulkmid.t, method="unifrac", weighted=F)
ordination = ordinate(psbulkmid.t, method="PCoA", distance=psbulkmidtunifrac_dist)
plot_ordination(psbulkmid.t, ordination, color = "Plot",) + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psbulkmid.t.tree.data <- data.frame(sample_data(psbulkmid.t))
adonis2(psbulkmidtunifrac_dist ~Plot, psbulkmid.t.tree.data)

#betadispersion
psbulkmidtu.bdisper <- betadisper(psbulkmidtunifrac_dist, psbulkmid.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulkmidtu.bdisper)
permutest(psbulkmidtu.bdisper, pairwise = TRUE)


#Weighted Bulk Mid
psbulkmidtunifrac_dist <- phyloseq::distance(psbulkmid.t, method="unifrac", weighted=T)
ordination = ordinate(psbulkmid.t, method="PCoA", distance=psbulkmidtunifrac_dist)
plot_ordination(psbulkmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psbulkmid.t.tree.data <- data.frame(sample_data(psbulkmid.t))
adonis2(psbulkmidtunifrac_dist ~Plot, psbulkmid.t.tree.data)

#betadispersion
psbulkmidtu.bdisper <- betadisper(psbulkmidtunifrac_dist, psbulkmid.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulkmidtu.bdisper)
permutest(psbulkmidtu.bdisper, pairwise = TRUE)


#Beta diversity- Unifrac Bulk Soil Fall season samples ----
psbulkfall.t.tree <- rtree(ntaxa(psbulkfall.t), rooted = TRUE, tip.label = taxa_names(psbulkfall.t))
#plot(psbulkfall.t.tree)
psbulkfall.t <- merge_phyloseq(psbulkfall.t, psbulkfall.t.tree)

#Unweighted Bulk Fall
psbulkfalltunifrac_dist <- phyloseq::distance(psbulkfall.t, method="unifrac", weighted=F)
ordination = ordinate(psbulkfall.t, method="PCoA", distance=psbulkfalltunifrac_dist)
plot_ordination(psbulkfall.t, ordination, color = "Plot",) + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psbulkfall.t.tree.data <- data.frame(sample_data(psbulkfall.t))
adonis2(psbulkfalltunifrac_dist ~Plot, psbulkfall.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psbulkfalltunifrac_dist, psbulkfall.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
psbulkfalltu.bdisper <- betadisper(psbulkfalltunifrac_dist, psbulkfall.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulkfalltu.bdisper)
permutest(psbulkfalltu.bdisper, pairwise = TRUE)



#Weighted Bulk Fall
psbulkfalltunifrac_dist <- phyloseq::distance(psbulkfall.t, method="unifrac", weighted=T)
ordination = ordinate(psbulkfall.t, method="PCoA", distance=psbulkfalltunifrac_dist)
plot_ordination(psbulkfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psbulkfall.t.tree.data <- data.frame(sample_data(psbulkfall.t))
adonis2(psbulkfalltunifrac_dist ~Plot, psbulkfall.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psbulkfalltunifrac_dist, psbulkfall.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
psbulkfalltu.bdisper <- betadisper(psbulkfalltunifrac_dist, psbulkfall.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulkfalltu.bdisper)
permutest(psbulkfalltu.bdisper, pairwise = TRUE)



#Beta diversity- Unifrac Cyst ----
pscyst.t.tree <- rtree(ntaxa(pscyst.t), rooted = TRUE, tip.label = taxa_names(pscyst.t))
#plot(pscyst.t.tree)
pscyst.t <- merge_phyloseq(pscyst.t, pscyst.t.tree)

#Unweighted
pscysttunifrac_dist <- phyloseq::distance(pscyst.t, method="unifrac", weighted=F)
ordination = ordinate(pscyst.t, method="PCoA", distance=pscysttunifrac_dist)
#plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pscyst.t.tree.data <- data.frame(sample_data(pscyst.t))
adonis2(pscysttunifrac_dist ~Plot*Season*Num_Cysts, pscyst.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pscysttunifrac_dist, pscyst.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
pscysttu.bdisper <- betadisper(pscysttunifrac_dist, pscyst.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscysttu.bdisper)
permutest(pscysttu.bdisper, pairwise = TRUE)


#Weighted
pscysttunifrac_dist <- phyloseq::distance(pscyst.t, method="unifrac", weighted=T)
ordination = ordinate(pscyst.t, method="PCoA", distance=pscysttunifrac_dist)
#plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pscyst.t.tree.data <- data.frame(sample_data(pscyst.t))
adonis2(pscysttunifrac_dist ~Plot*Season*Num_Cysts, pscyst.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pscysttunifrac_dist, pscyst.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
pscysttu.bdisper <- betadisper(pscysttunifrac_dist, pscyst.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscysttu.bdisper)
permutest(pscysttu.bdisper, pairwise = TRUE)

#Beta diversity- Unifrac Cysts mid season samples ----
pscystmid.t.tree <- rtree(ntaxa(pscystmid.t), rooted = TRUE, tip.label = taxa_names(pscystmid.t))
#plot(pscystmid.t.tree)
pscystmid.t <- merge_phyloseq(pscystmid.t, pscystmid.t.tree)

#Unweighted Cysts mid
pscystmidtunifrac_dist <- phyloseq::distance(pscystmid.t, method="unifrac", weighted=F)
ordination = ordinate(pscystmid.t, method="PCoA", distance=pscystmidtunifrac_dist)
#plot_ordination(pscystmid.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pscystmid.t.tree.data <- data.frame(sample_data(pscystmid.t))
adonis2(pscystmidtunifrac_dist ~Plot*Num_Cysts, pscystmid.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(pscystmidtunifrac_dist, pscystmid.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
pscystmidtu.bdisper <- betadisper(pscystmidtunifrac_dist, pscystmid.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscystmidtu.bdisper)
permutest(pscystmidtu.bdisper, pairwise = TRUE)


#Weighted Cysts mid
pscystmidtunifrac_dist <- phyloseq::distance(pscystmid.t, method="unifrac", weighted=T)
ordination = ordinate(pscystmid.t, method="PCoA", distance=pscystmidtunifrac_dist)
#plot_ordination(pscystmid.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pscystmid.t.tree.data <- data.frame(sample_data(pscystmid.t))
adonis2(pscystmidtunifrac_dist ~Plot*Num_Cysts, pscystmid.t.tree.data)

#betadispersion
pscystmidtu.bdisper <- betadisper(pscystmidtunifrac_dist, pscystmid.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscystmidtu.bdisper)
permutest(pscystmidtu.bdisper, pairwise = TRUE)


#Beta diversity- Unifrac Cysts fall season samples ----
pscystfall.t.tree <- rtree(ntaxa(pscystfall.t), rooted = TRUE, tip.label = taxa_names(pscystfall.t))
#plot(pscystfall.t.tree)
pscystfall.t <- merge_phyloseq(pscystfall.t, pscystfall.t.tree)

#Unweighted Cysts Fall
pscystfalltunifrac_dist <- phyloseq::distance(pscystfall.t, method="unifrac", weighted=F)
ordination = ordinate(pscystfall.t, method="PCoA", distance=pscystfalltunifrac_dist)
#plot_ordination(pscystfall.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pscystfall.t.tree.data <- data.frame(sample_data(pscystfall.t))
adonis2(pscystfalltunifrac_dist ~Plot, pscystfall.t.tree.data)

#betadispersion
pscystfalltu.bdisper <- betadisper(pscystfalltunifrac_dist, pscystfall.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscystfalltu.bdisper)
permutest(pscystfalltu.bdisper, pairwise = TRUE)


#Weighted Cysts Fall
pscystfalltunifrac_dist <- phyloseq::distance(pscystfall.t, method="unifrac", weighted=T)
ordination = ordinate(pscystfall.t, method="PCoA", distance=pscystfalltunifrac_dist)
#plot_ordination(pscystfall.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
pscystfall.t.tree.data <- data.frame(sample_data(pscystfall.t))
adonis2(pscystfalltunifrac_dist ~Plot, pscystfall.t.tree.data)

#betadispersion
pscystfalltu.bdisper <- betadisper(pscystfalltunifrac_dist, pscystfall.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscystfalltu.bdisper)
permutest(pscystfalltu.bdisper, pairwise = TRUE)


#Beta diversity- Unifrac Roots ----
psroot.t.tree <- rtree(ntaxa(psroot.t), rooted = TRUE, tip.label = taxa_names(psroot.t))
#plot(psroot.t.tree)
psroot.t <- merge_phyloseq(psroot.t, psroot.t.tree)

#Unweighted Root
psroottunifrac_dist <- phyloseq::distance(psroot.t, method="unifrac", weighted=F)
ordination = ordinate(psroot.t, method="PCoA", distance=psroottunifrac_dist)
plot_ordination(psroot.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psroot.t.tree.data <- data.frame(sample_data(psroot.t))
adonis2(psroottunifrac_dist ~Plot, psroot.t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psroottunifrac_dist, psroot.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
psroottu.bdisper <- betadisper(psroottunifrac_dist, psroot.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psroottu.bdisper)
permutest(psroottu.bdisper, pairwise = TRUE)


#Weighted Root
psroottunifrac_dist <- phyloseq::distance(psroot.t, method="unifrac", weighted=T)
ordination = ordinate(psroot.t, method="PCoA", distance=psroottunifrac_dist)
plot_ordination(psroot.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)

#statistics
#permanova
psroot.t.tree.data <- data.frame(sample_data(psroot.t))
adonis2(psroottunifrac_dist ~Plot, psroot.t.tree.data)

#pairwise posthoc
pairwiseAdonis::pairwise.adonis(psroottunifrac_dist, psroot.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)


#betadispersion
psroottu.bdisper <- betadisper(psroottunifrac_dist, psroot.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psroottu.bdisper)
permutest(psroottu.bdisper, pairwise = TRUE)


#Taxonomy and Prevalence Heatmaps ----
# Split 16S bulk soil into individual objects ----
psbulkSA.16S <- subset_samples(psbulk, Plot == "Sa")
psbulkS3.16S <- subset_samples(psbulk, Plot == "S3")
psbulkSS.16S <- subset_samples(psbulk, Plot == "Ss")
psbulkSR.16S <- subset_samples(psbulk, Plot == "Sr")

#split into mid season
psbulkSAmid.16S <- subset_samples(psbulkSA.16S, Season == "Mid")
psbulkS3mid.16S <- subset_samples(psbulkS3.16S, Season == "Mid")
psbulkSSmid.16S <- subset_samples(psbulkSS.16S, Season == "Mid")
psbulkSRmid.16S <- subset_samples(psbulkSR.16S, Season == "Mid")

#split into fall season
psbulkSAfall.16S <- subset_samples(psbulkSA.16S, Season == "Fall")
psbulkS3fall.16S <- subset_samples(psbulkS3.16S, Season == "Fall")
psbulkSSfall.16S <- subset_samples(psbulkSS.16S, Season == "Fall")
psbulkSRfall.16S <- subset_samples(psbulkSR.16S, Season == "Fall")

#remove taxa that have 0 sums.
psbulkSAmid.16S <- prune_taxa(taxa_sums(psbulkSAmid.16S) > 0, psbulkSAmid.16S)
psbulkS3mid.16S <- prune_taxa(taxa_sums(psbulkS3mid.16S) > 0, psbulkS3mid.16S)
psbulkSRmid.16S <- prune_taxa(taxa_sums(psbulkSRmid.16S) > 0, psbulkSRmid.16S)
psbulkSSmid.16S <- prune_taxa(taxa_sums(psbulkSSmid.16S) > 0, psbulkSSmid.16S)

psbulkSAfall.16S <- prune_taxa(taxa_sums(psbulkSAfall.16S) > 0, psbulkSAfall.16S)
psbulkS3fall.16S <- prune_taxa(taxa_sums(psbulkS3fall.16S) > 0, psbulkS3fall.16S)
psbulkSRfall.16S <- prune_taxa(taxa_sums(psbulkSRfall.16S) > 0, psbulkSRfall.16S)
psbulkSSfall.16S <- prune_taxa(taxa_sums(psbulkSSfall.16S) > 0, psbulkSSfall.16S)


# Split 16S cysts into individual objects ----
pscystSA.16S <- subset_samples(pscyst, Plot == "Sa")
pscystS3.16S <- subset_samples(pscyst, Plot == "S3")
pscystSS.16S <- subset_samples(pscyst, Plot == "Ss")
pscystSR.16S <- subset_samples(pscyst, Plot == "Sr")

#split into mid season
pscystSAmid.16S <- subset_samples(pscystSA.16S, Season == "Mid")
pscystS3mid.16S <- subset_samples(pscystS3.16S, Season == "Mid")
pscystSSmid.16S <- subset_samples(pscystSS.16S, Season == "Mid")
pscystSRmid.16S <- subset_samples(pscystSR.16S, Season == "Mid")
#split into fall season
pscystSAfall.16S <- subset_samples(pscystSA.16S, Season == "Fall")
pscystS3fall.16S <- subset_samples(pscystS3.16S, Season == "Fall")
pscystSSfall.16S <- subset_samples(pscystSS.16S, Season == "Fall")
pscystSRfall.16S <- subset_samples(pscystSR.16S, Season == "Fall")

#remove taxa that have 0 sums.
pscystSAmid.16S <- prune_taxa(taxa_sums(pscystSAmid.16S) > 0, pscystSAmid.16S)
pscystS3mid.16S <- prune_taxa(taxa_sums(pscystS3mid.16S) > 0, pscystS3mid.16S)
pscystSRmid.16S <- prune_taxa(taxa_sums(pscystSRmid.16S) > 0, pscystSRmid.16S)
pscystSSmid.16S <- prune_taxa(taxa_sums(pscystSSmid.16S) > 0, pscystSSmid.16S)

pscystSAfall.16S <- prune_taxa(taxa_sums(pscystSAfall.16S) > 0, pscystSAfall.16S)
pscystS3fall.16S <- prune_taxa(taxa_sums(pscystS3fall.16S) > 0, pscystS3fall.16S)
pscystSRfall.16S <- prune_taxa(taxa_sums(pscystSRfall.16S) > 0, pscystSRfall.16S)
pscystSSfall.16S <- prune_taxa(taxa_sums(pscystSSfall.16S) > 0, pscystSSfall.16S)

# Split 16S Rhizosphere into individual objects ----
psrhizSA.16S <- subset_samples(psrhiz, Plot == "Sa")
psrhizS3.16S <- subset_samples(psrhiz, Plot == "S3")
psrhizSS.16S <- subset_samples(psrhiz, Plot == "Ss")
psrhizSR.16S <- subset_samples(psrhiz, Plot == "Sr")

#remove taxa that have 0 sums.
psrhizSA.16S <- prune_taxa(taxa_sums(psrhizSA.16S) > 0, psrhizSA.16S)
psrhizS3.16S <- prune_taxa(taxa_sums(psrhizS3.16S) > 0, psrhizS3.16S)
psrhizSR.16S <- prune_taxa(taxa_sums(psrhizSR.16S) > 0, psrhizSR.16S)
psrhizSS.16S <- prune_taxa(taxa_sums(psrhizSS.16S) > 0, psrhizSS.16S)

# Split 16S Root into individual objects ----
psrootSA.16S <- subset_samples(psroot, Plot == "Sa")
psrootS3.16S <- subset_samples(psroot, Plot == "S3")
psrootSS.16S <- subset_samples(psroot, Plot == "Ss")
psrootSR.16S <- subset_samples(psroot, Plot == "Sr")

#remove taxa that have 0 sums.
psrootSA.16S <- prune_taxa(taxa_sums(psrootSA.16S) > 0, psrootSA.16S)
psrootS3.16S <- prune_taxa(taxa_sums(psrootS3.16S) > 0, psrootS3.16S)
psrootSR.16S <- prune_taxa(taxa_sums(psrootSR.16S) > 0, psrootSR.16S)
psrootSS.16S <- prune_taxa(taxa_sums(psrootSS.16S) > 0, psrootSS.16S)



# Cysts 16S plots- taxa with 2% abundance or higher ----

#To make taxa be listed from high abundance to low abundance, ggplot command instead of fill = Genus
#replace with: factor(Genus, levels = c(setdiff(Genus, "Other"), "Other"))

#S3 mid season
pscystS3mid.16S.rab <- transform_sample_counts(pscystS3mid.16S, function(x) x/sum(x))
filterpscystS3mid.16S.rab <- phyloseq::genefilter_sample(pscystS3mid.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3mid.16Score.rab <- prune_taxa(filterpscystS3mid.16S.rab, pscystS3mid.16S.rab)
cystS3mid.16Score.genus.rab <- tax_glom(cystS3mid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#four NA are being combined. Need to make sure they stay separate
cystS3mid.16Score.genus.rab <-  renameTaxa(cystS3mid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystS3mid.16Score.genus.rab.df <- psmelt(cystS3mid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3mid.16Score.genus.rab.df$Genus <- as.character(cystS3mid.16Score.genus.rab.df$Genus)
ggplot(cystS3mid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA mid season
pscystSAmid.16S.rab <- transform_sample_counts(pscystSAmid.16S, function(x) x/sum(x))
filterpscystSAmid.16S.rab <- phyloseq::genefilter_sample(pscystSAmid.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSAmid.16Score.rab <- prune_taxa(filterpscystSAmid.16S.rab, pscystSAmid.16S.rab)
cystSAmid.16Score.genus.rab <- tax_glom(cystSAmid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSAmid.16Score.genus.rab.df <- psmelt(cystSAmid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSAmid.16Score.genus.rab.df$Genus <- as.character(cystSAmid.16Score.genus.rab.df$Genus)
ggplot(cystSAmid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR mid season
pscystSRmid.16S.rab <- transform_sample_counts(pscystSRmid.16S, function(x) x/sum(x))
filterpscystSRmid.16S.rab <- phyloseq::genefilter_sample(pscystSRmid.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSRmid.16Score.rab <- prune_taxa(filterpscystSRmid.16S.rab, pscystSRmid.16S.rab)
cystSRmid.16Score.genus.rab <- tax_glom(cystSRmid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#two NA are being combined. Need to make sure they stay separate
cystSRmid.16Score.genus.rab <-  renameTaxa(cystSRmid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSRmid.16Score.genus.rab.df <- psmelt(cystSRmid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSRmid.16Score.genus.rab.df$Genus <- as.character(cystSRmid.16Score.genus.rab.df$Genus)
ggplot(cystSRmid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS mid season
pscystSSmid.16S.rab <- transform_sample_counts(pscystSSmid.16S, function(x) x/sum(x))
filterpscystSSmid.16S.rab <- phyloseq::genefilter_sample(pscystSSmid.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSSmid.16Score.rab <- prune_taxa(filterpscystSSmid.16S.rab, pscystSSmid.16S.rab)
cystSSmid.16Score.genus.rab <- tax_glom(cystSSmid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#two NA are being combined. Need to make sure they stay separate
cystSSmid.16Score.genus.rab <-  renameTaxa(cystSSmid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSSmid.16Score.genus.rab.df <- psmelt(cystSSmid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSSmid.16Score.genus.rab.df$Genus <- as.character(cystSSmid.16Score.genus.rab.df$Genus)
ggplot(cystSSmid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#S3 fall season
pscystS3fall.16S.rab <- transform_sample_counts(pscystS3fall.16S, function(x) x/sum(x))
filterpscystS3fall.16S.rab <- phyloseq::genefilter_sample(pscystS3fall.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3fall.16Score.rab <- prune_taxa(filterpscystS3fall.16S.rab, pscystS3fall.16S.rab)
cystS3fall.16Score.genus.rab <- tax_glom(cystS3fall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#one NA was unlabeled.
cystS3fall.16Score.genus.rab <-  renameTaxa(cystS3fall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystS3fall.16Score.genus.rab.df <- psmelt(cystS3fall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3fall.16Score.genus.rab.df$Genus <- as.character(cystS3fall.16Score.genus.rab.df$Genus)
ggplot(cystS3fall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA fall season
pscystSAfall.16S.rab <- transform_sample_counts(pscystSAfall.16S, function(x) x/sum(x))
filterpscystSAfall.16S.rab <- phyloseq::genefilter_sample(pscystSAfall.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSAfall.16Score.rab <- prune_taxa(filterpscystSAfall.16S.rab, pscystSAfall.16S.rab)
cystSAfall.16Score.genus.rab <- tax_glom(cystSAfall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#one NA was unlabeled.
cystSAfall.16Score.genus.rab <-  renameTaxa(cystSAfall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSAfall.16Score.genus.rab.df <- psmelt(cystSAfall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSAfall.16Score.genus.rab.df$Genus <- as.character(cystSAfall.16Score.genus.rab.df$Genus)
ggplot(cystSAfall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR fall season
pscystSRfall.16S.rab <- transform_sample_counts(pscystSRfall.16S, function(x) x/sum(x))
filterpscystSRfall.16S.rab <- phyloseq::genefilter_sample(pscystSRfall.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSRfall.16Score.rab <- prune_taxa(filterpscystSRfall.16S.rab, pscystSRfall.16S.rab)
cystSRfall.16Score.genus.rab <- tax_glom(cystSRfall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#two NA are being combined. Need to make sure they stay separate
cystSRfall.16Score.genus.rab <-  renameTaxa(cystSRfall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSRfall.16Score.genus.rab.df <- psmelt(cystSRfall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSRfall.16Score.genus.rab.df$Genus <- as.character(cystSRfall.16Score.genus.rab.df$Genus)
ggplot(cystSRfall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS fall season
pscystSSfall.16S.rab <- transform_sample_counts(pscystSSfall.16S, function(x) x/sum(x))
filterpscystSSfall.16S.rab <- phyloseq::genefilter_sample(pscystSSfall.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSSfall.16Score.rab <- prune_taxa(filterpscystSSfall.16S.rab, pscystSSfall.16S.rab)
cystSSfall.16Score.genus.rab <- tax_glom(cystSSfall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#four NA are being combined. Need to make sure they stay separate
cystSSfall.16Score.genus.rab <-  renameTaxa(cystSSfall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSSfall.16Score.genus.rab.df <- psmelt(cystSSfall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSSfall.16Score.genus.rab.df$Genus <- as.character(cystSSfall.16Score.genus.rab.df$Genus)
ggplot(cystSSfall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

# Merge the Cyst 2% abundance or higher dataframes together then plot heatmap. Seasons are split----

allcystmid.16Score.genus.df <- rbind(cystS3mid.16Score.genus.rab.df, cystSAmid.16Score.genus.rab.df, cystSRmid.16Score.genus.rab.df, cystSSmid.16Score.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allcystmid.16Score.genus.df, "allcystmid.16Score.genus.csv")
#import back in
allcystmid.16Score.genus.df <- read.csv("allcystmid.16Score.genus.csv")

#ggplot(allcystmid.16Score.genus.df, aes(x = Sample, y = Genus, fill = Abundance)) + geom_tile(color = "black", 
#) + theme(axis.text.x = element_blank()) + facet_grid(~Plot, scales = "free", space = "free") + scale_fill_gradient(low = "yellow", high = "blue")
#ggplot(allcystmid.16Score.genus.df, aes(x = Sample, y = Genus, fill = Abundance)) + geom_tile(color = "black"
#) + theme(axis.text.x = element_blank()) + facet_grid(~Plot, scales = "free", space = "free") + scale_fill_gradientn(colors = hcl.colors(10, "Oslo"))

#make legend max, 53%
ggplot(allcystmid.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = 0.27
                                   , limits=c(0.0, 0.53), breaks=seq(0.0, 0.53, by=0.10))



allcystfall.16Score.genus.df <- rbind(cystS3fall.16Score.genus.rab.df, cystSAfall.16Score.genus.rab.df, cystSRfall.16Score.genus.rab.df, cystSSfall.16Score.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allcystfall.16Score.genus.df, "allcystfall.16Score.genus.csv")
#import back in
allcystfall.16Score.genus.df <- read.csv("allcystfall.16Score.genus.csv")

#max legend 54%
ggplot(allcystfall.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black"
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"),axis.text.y.left = element_text(size =11),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = 0.27, limits=c(0.0, 0.54), breaks=seq(0.0, 0.54 ,by=0.09))


#ggplot(allcystfall.16Score.genus.df, aes(x = Sample, y = Genus, fill = Abundance)) + geom_tile(color = "black", ) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#axis.text.y.left = element_text(size =11)) + facet_grid(~Plot, scales = "free", space = "free") + scale_fill_gradient(low = "blue",  high = "orange")




# Cysts 16S- treatments- not split by season- 2% abundance or higher----

#pscystSA.16S <- subset_samples(pscyst, Plot == "Sa")
#pscystS3.16S <- subset_samples(pscyst, Plot == "S3")
#pscystSS.16S <- subset_samples(pscyst, Plot == "Ss")
#pscystSR.16S <- subset_samples(pscyst, Plot == "Sr")

#S3  
pscystS3.16S.rab <- transform_sample_counts(pscystS3.16S, function(x) x/sum(x))
filterpscystS3.16S.rab <- phyloseq::genefilter_sample(pscystS3.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3.16Score.rab <- prune_taxa(filterpscystS3.16S.rab, pscystS3.16S.rab)
cystS3.16Score.genus.rab <- tax_glom(cystS3.16Score.rab, taxrank = "Genus", NArm = FALSE)
#4 NA are being combined. Need to make sure they stay separate
cystS3.16Score.genus.rab <-  renameTaxa(cystS3.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystS3.16Score.genus.rab.df <- psmelt(cystS3.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3.16Score.genus.rab.df$Genus <- as.character(cystS3.16Score.genus.rab.df$Genus)
ggplot(cystS3.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA  
pscystSA.16S.rab <- transform_sample_counts(pscystSA.16S, function(x) x/sum(x))
filterpscystSA.16S.rab <- phyloseq::genefilter_sample(pscystSA.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSA.16Score.rab <- prune_taxa(filterpscystSA.16S.rab, pscystSA.16S.rab)
cystSA.16Score.genus.rab <- tax_glom(cystSA.16Score.rab, taxrank = "Genus", NArm = FALSE)
#1 NA was being combined. Need to make sure they stay separate
cystSA.16Score.genus.rab <-  renameTaxa(cystSA.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSA.16Score.genus.rab.df <- psmelt(cystSA.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSA.16Score.genus.rab.df$Genus <- as.character(cystSA.16Score.genus.rab.df$Genus)
ggplot(cystSA.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR  
pscystSR.16S.rab <- transform_sample_counts(pscystSR.16S, function(x) x/sum(x))
filterpscystSR.16S.rab <- phyloseq::genefilter_sample(pscystSR.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSR.16Score.rab <- prune_taxa(filterpscystSR.16S.rab, pscystSR.16S.rab)
cystSR.16Score.genus.rab <- tax_glom(cystSR.16Score.rab, taxrank = "Genus", NArm = FALSE)
#3 NA are being combined. Need to make sure they stay separate
cystSR.16Score.genus.rab <-  renameTaxa(cystSR.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSR.16Score.genus.rab.df <- psmelt(cystSR.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSR.16Score.genus.rab.df$Genus <- as.character(cystSR.16Score.genus.rab.df$Genus)
ggplot(cystSR.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS  
pscystSS.16S.rab <- transform_sample_counts(pscystSS.16S, function(x) x/sum(x))
filterpscystSS.16S.rab <- phyloseq::genefilter_sample(pscystSS.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSS.16Score.rab <- prune_taxa(filterpscystSS.16S.rab, pscystSS.16S.rab)
cystSS.16Score.genus.rab <- tax_glom(cystSS.16Score.rab, taxrank = "Genus", NArm = FALSE)
#4 NA are being combined. Need to make sure they stay separate
cystSS.16Score.genus.rab <-  renameTaxa(cystSS.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
cystSS.16Score.genus.rab.df <- psmelt(cystSS.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
cystSS.16Score.genus.rab.df$Genus <- as.character(cystSS.16Score.genus.rab.df$Genus)
ggplot(cystSS.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")



# Merge Cyst 2% abundance or higher dataframes together- plot heatmap- not split by season. ----
allcysts.16Score.genus.df <- rbind(cystS3.16Score.genus.rab.df, cystSA.16Score.genus.rab.df, cystSR.16Score.genus.rab.df, cystSS.16Score.genus.rab.df)

#fix names of family taxa to have the same name- I removed the numbers introduced by the renaming. 
#write.csv(allcysts.16Score.genus.df, "allcysts.16Score.genus.csv")

#import back in
allcysts.16Score.genus.df <- read.csv("allcysts.16Score.genus.csv")

#make mid season first
allcysts.16Score.genus.df$Season <- factor(allcysts.16Score.genus.df$Season, levels = c("Mid", "Fall"))


#added prevlence and total abundance to the same file, just renamed it in case it fails.
allcysts.16Score.genus.prev.df <- read.csv("allcysts.16Score.genus.prev.csv")

#reorder seasons
allcysts.16Score.genus.prev.df$Season <- factor(allcysts.16Score.genus.prev.df$Season, levels = c("Mid", "Fall"))

#need to change the values depending on how you want to view it. 
ggplot(allcysts.16Score.genus.prev.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
           , breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.537), values = c("0", "0.2", "1"))

#break


# Bulk 16S- treatments- not split by season - 1% or higher----

#psbulkSA.16S <- subset_samples(psbulk, Plot == "Sa")
#psbulkS3.16S <- subset_samples(psbulk, Plot == "S3")
#psbulkSS.16S <- subset_samples(psbulk, Plot == "Ss")
#psbulkSR.16S <- subset_samples(psbulk, Plot == "Sr")


#S3  
psbulkS3.16S.rab <- transform_sample_counts(psbulkS3.16S, function(x) x/sum(x))
filterpsbulkS3.16S.rab <- phyloseq::genefilter_sample(psbulkS3.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkS3.16Score.rab <- prune_taxa(filterpsbulkS3.16S.rab, psbulkS3.16S.rab)
bulkS3.16Score.genus.rab <- tax_glom(bulkS3.16Score.rab, taxrank = "Genus", NArm = FALSE)
#3 NA are being combined. Need to make sure they stay separate
bulkS3.16Score.genus.rab <-  renameTaxa(bulkS3.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkS3.16Score.genus.rab.df <- psmelt(bulkS3.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3.16Score.genus.rab.df$Genus <- as.character(bulkS3.16Score.genus.rab.df$Genus)
ggplot(bulkS3.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA  
psbulkSA.16S.rab <- transform_sample_counts(psbulkSA.16S, function(x) x/sum(x))
filterpsbulkSA.16S.rab <- phyloseq::genefilter_sample(psbulkSA.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSA.16Score.rab <- prune_taxa(filterpsbulkSA.16S.rab, psbulkSA.16S.rab)
bulkSA.16Score.genus.rab <- tax_glom(bulkSA.16Score.rab, taxrank = "Genus", NArm = FALSE)
#6 NA was being combined. Need to make sure they stay separate
bulkSA.16Score.genus.rab <-  renameTaxa(bulkSA.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSA.16Score.genus.rab.df <- psmelt(bulkSA.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSA.16Score.genus.rab.df$Genus <- as.character(bulkSA.16Score.genus.rab.df$Genus)
ggplot(bulkSA.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR  
psbulkSR.16S.rab <- transform_sample_counts(psbulkSR.16S, function(x) x/sum(x))
filterpsbulkSR.16S.rab <- phyloseq::genefilter_sample(psbulkSR.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSR.16Score.rab <- prune_taxa(filterpsbulkSR.16S.rab, psbulkSR.16S.rab)
bulkSR.16Score.genus.rab <- tax_glom(bulkSR.16Score.rab, taxrank = "Genus", NArm = FALSE)
#4 NA are being combined. Need to make sure they stay separate
bulkSR.16Score.genus.rab <-  renameTaxa(bulkSR.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSR.16Score.genus.rab.df <- psmelt(bulkSR.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSR.16Score.genus.rab.df$Genus <- as.character(bulkSR.16Score.genus.rab.df$Genus)
ggplot(bulkSR.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS  
psbulkSS.16S.rab <- transform_sample_counts(psbulkSS.16S, function(x) x/sum(x))
filterpsbulkSS.16S.rab <- phyloseq::genefilter_sample(psbulkSS.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSS.16Score.rab <- prune_taxa(filterpsbulkSS.16S.rab, psbulkSS.16S.rab)
bulkSS.16Score.genus.rab <- tax_glom(bulkSS.16Score.rab, taxrank = "Genus", NArm = FALSE)
#4 NA are being combined. Need to make sure they stay separate
bulkSS.16Score.genus.rab <-  renameTaxa(bulkSS.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSS.16Score.genus.rab.df <- psmelt(bulkSS.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSS.16Score.genus.rab.df$Genus <- as.character(bulkSS.16Score.genus.rab.df$Genus)
ggplot(bulkSS.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

#break


# Bulk 16S- Heatmap of 1% or higher- not split by season ----

allbulk.16Score.genus.df <- rbind(bulkS3.16Score.genus.rab.df, bulkSA.16Score.genus.rab.df, bulkSR.16Score.genus.rab.df, bulkSS.16Score.genus.rab.df)

write.csv(allbulk.16Score.genus.df, "allbulk.16Score.genus.csv")
#fix names of family taxa to have the same name.
#added prevlence and total abundance too
#import back in
allbulk.16Score.genus.df <- read.csv("fix.allbulk.16Score.genus.csv")


#make mid season first
allbulk.16Score.genus.df$Season <- factor(allbulk.16Score.genus.df$Season, levels = c("Mid", "Fall"))


#WINNER
#max % is 29.8
#need to change the values depending on how you want to view it. Values are 0 to 1.  
ggplot(allbulk.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.1, 0.2, 0.298), values = c("0", "0.35", "1"))


#break

# Bulk Soil 16S plots- split by season- taxa with 1% abundance or higher ----

#To make taxa be listed from high abundance to low abundance
#in ggplot command instead of fill = Genus
#replace with: factor(Genus, levels = c(setdiff(Genus, "Other"), "Other"))

#S3 mid season
psbulkS3mid.16S.rab <- transform_sample_counts(psbulkS3mid.16S, function(x) x/sum(x))
filterpsbulkS3mid.16S.rab <- phyloseq::genefilter_sample(psbulkS3mid.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkS3mid.16Score.rab <- prune_taxa(filterpsbulkS3mid.16S.rab, psbulkS3mid.16S.rab)
bulkS3mid.16Score.genus.rab <- tax_glom(bulkS3mid.16Score.rab, taxrank = "Genus", NArm = FALSE)
# NA are being combined. Need to make sure they stay separate
bulkS3mid.16Score.genus.rab <-  renameTaxa(bulkS3mid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkS3mid.16Score.genus.rab.df <- psmelt(bulkS3mid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3mid.16Score.genus.rab.df$Genus <- as.character(bulkS3mid.16Score.genus.rab.df$Genus)
ggplot(bulkS3mid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA mid season
psbulkSAmid.16S.rab <- transform_sample_counts(psbulkSAmid.16S, function(x) x/sum(x))
filterpsbulkSAmid.16S.rab <- phyloseq::genefilter_sample(psbulkSAmid.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSAmid.16Score.rab <- prune_taxa(filterpsbulkSAmid.16S.rab, psbulkSAmid.16S.rab)
bulkSAmid.16Score.genus.rab <- tax_glom(bulkSAmid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#2 NAs are being combined. Need to make sure they stay separate
bulkSAmid.16Score.genus.rab <-  renameTaxa(bulkSAmid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe
bulkSAmid.16Score.genus.rab.df <- psmelt(bulkSAmid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSAmid.16Score.genus.rab.df$Genus <- as.character(bulkSAmid.16Score.genus.rab.df$Genus)
ggplot(bulkSAmid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR mid season
psbulkSRmid.16S.rab <- transform_sample_counts(psbulkSRmid.16S, function(x) x/sum(x))
filterpsbulkSRmid.16S.rab <- phyloseq::genefilter_sample(psbulkSRmid.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSRmid.16Score.rab <- prune_taxa(filterpsbulkSRmid.16S.rab, psbulkSRmid.16S.rab)
bulkSRmid.16Score.genus.rab <- tax_glom(bulkSRmid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#3 NA are being combined. Need to make sure they stay separate
bulkSRmid.16Score.genus.rab <-  renameTaxa(bulkSRmid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSRmid.16Score.genus.rab.df <- psmelt(bulkSRmid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSRmid.16Score.genus.rab.df$Genus <- as.character(bulkSRmid.16Score.genus.rab.df$Genus)
ggplot(bulkSRmid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS mid season
psbulkSSmid.16S.rab <- transform_sample_counts(psbulkSSmid.16S, function(x) x/sum(x))
filterpsbulkSSmid.16S.rab <- phyloseq::genefilter_sample(psbulkSSmid.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSSmid.16Score.rab <- prune_taxa(filterpsbulkSSmid.16S.rab, psbulkSSmid.16S.rab)
bulkSSmid.16Score.genus.rab <- tax_glom(bulkSSmid.16Score.rab, taxrank = "Genus", NArm = FALSE)
#3 NA are being combined. Need to make sure they stay separate
bulkSSmid.16Score.genus.rab <-  renameTaxa(bulkSSmid.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSSmid.16Score.genus.rab.df <- psmelt(bulkSSmid.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSSmid.16Score.genus.rab.df$Genus <- as.character(bulkSSmid.16Score.genus.rab.df$Genus)
ggplot(bulkSSmid.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#S3 fall season
psbulkS3fall.16S.rab <- transform_sample_counts(psbulkS3fall.16S, function(x) x/sum(x))
filterpsbulkS3fall.16S.rab <- phyloseq::genefilter_sample(psbulkS3fall.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkS3fall.16Score.rab <- prune_taxa(filterpsbulkS3fall.16S.rab, psbulkS3fall.16S.rab)
bulkS3fall.16Score.genus.rab <- tax_glom(bulkS3fall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#3 NA was unlabeled.
bulkS3fall.16Score.genus.rab <-  renameTaxa(bulkS3fall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkS3fall.16Score.genus.rab.df <- psmelt(bulkS3fall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3fall.16Score.genus.rab.df$Genus <- as.character(bulkS3fall.16Score.genus.rab.df$Genus)
ggplot(bulkS3fall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA fall season
psbulkSAfall.16S.rab <- transform_sample_counts(psbulkSAfall.16S, function(x) x/sum(x))
filterpsbulkSAfall.16S.rab <- phyloseq::genefilter_sample(psbulkSAfall.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSAfall.16Score.rab <- prune_taxa(filterpsbulkSAfall.16S.rab, psbulkSAfall.16S.rab)
bulkSAfall.16Score.genus.rab <- tax_glom(bulkSAfall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#6 NA was unlabeled.
bulkSAfall.16Score.genus.rab <-  renameTaxa(bulkSAfall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSAfall.16Score.genus.rab.df <- psmelt(bulkSAfall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSAfall.16Score.genus.rab.df$Genus <- as.character(bulkSAfall.16Score.genus.rab.df$Genus)
ggplot(bulkSAfall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR fall season
psbulkSRfall.16S.rab <- transform_sample_counts(psbulkSRfall.16S, function(x) x/sum(x))
filterpsbulkSRfall.16S.rab <- phyloseq::genefilter_sample(psbulkSRfall.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSRfall.16Score.rab <- prune_taxa(filterpsbulkSRfall.16S.rab, psbulkSRfall.16S.rab)
bulkSRfall.16Score.genus.rab <- tax_glom(bulkSRfall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#four NA are being combined. Need to make sure they stay separate
bulkSRfall.16Score.genus.rab <-  renameTaxa(bulkSRfall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSRfall.16Score.genus.rab.df <- psmelt(bulkSRfall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSRfall.16Score.genus.rab.df$Genus <- as.character(bulkSRfall.16Score.genus.rab.df$Genus)
ggplot(bulkSRfall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS fall season
psbulkSSfall.16S.rab <- transform_sample_counts(psbulkSSfall.16S, function(x) x/sum(x))
filterpsbulkSSfall.16S.rab <- phyloseq::genefilter_sample(psbulkSSfall.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSSfall.16Score.rab <- prune_taxa(filterpsbulkSSfall.16S.rab, psbulkSSfall.16S.rab)
bulkSSfall.16Score.genus.rab <- tax_glom(bulkSSfall.16Score.rab, taxrank = "Genus", NArm = FALSE)
#four NA are being combined. Need to make sure they stay separate
bulkSSfall.16Score.genus.rab <-  renameTaxa(bulkSSfall.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
bulkSSfall.16Score.genus.rab.df <- psmelt(bulkSSfall.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSSfall.16Score.genus.rab.df$Genus <- as.character(bulkSSfall.16Score.genus.rab.df$Genus)
ggplot(bulkSSfall.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

# Merge Bulk soil 16S dataframes together and make a heatmap of the 1% abundant or higher taxa ----

allbulkmid.16Score.genus.df <- rbind(bulkS3mid.16Score.genus.rab.df, bulkSAmid.16Score.genus.rab.df, bulkSRmid.16Score.genus.rab.df, bulkSSmid.16Score.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allbulkmid.16Score.genus.df, "allbulkmid.16Score.genus.csv")
#import back in
allbulkmid.16Score.genus.df <- read.csv("allbulkmid.16Score.genus.csv")


#label max legend as 20%
ggplot(allbulkmid.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"),axis.text.y.left = element_text(size =11),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = .1, limits=c(0.0, 0.20), breaks=seq(0.0, 0.20 ,by=0.05))



allbulkfall.16Score.genus.df <- rbind(bulkS3fall.16Score.genus.rab.df, bulkSAfall.16Score.genus.rab.df, bulkSRfall.16Score.genus.rab.df, bulkSSfall.16Score.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allbulkfall.16Score.genus.df, "allbulkfall.16Score.genus.csv")
#import back in
allbulkfall.16Score.genus.df <- read.csv("allbulkfall.16Score.genus.csv")

#label max percentage as 30%
ggplot(allbulkfall.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"),axis.text.y.left = element_text(size =11), 
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = 0.15, limits=c(0.0, 0.3), breaks=seq(0.0, 0.3 ,by=0.05))

#break


# Rhiz Soil 16S plots- taxa with 1% abundance or higher ----

#To make taxa be listed from high abundance to low abundance
#in ggplot command instead of fill = Genus
#replace with: factor(Genus, levels = c(setdiff(Genus, "Other"), "Other"))

#S3  season
psrhizS3.16S.rab <- transform_sample_counts(psrhizS3.16S, function(x) x/sum(x))
filterpsrhizS3.16S.rab <- phyloseq::genefilter_sample(psrhizS3.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizS3.16Score.rab <- prune_taxa(filterpsrhizS3.16S.rab, psrhizS3.16S.rab)
rhizS3.16Score.genus.rab <- tax_glom(rhizS3.16Score.rab, taxrank = "Genus", NArm = FALSE)
#three NA are being combined. Need to make sure they stay separate
rhizS3.16Score.genus.rab <-  renameTaxa(rhizS3.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rhizS3.16Score.genus.rab.df <- psmelt(rhizS3.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rhizS3.16Score.genus.rab.df$Genus <- as.character(rhizS3.16Score.genus.rab.df$Genus)
ggplot(rhizS3.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA  season
psrhizSA.16S.rab <- transform_sample_counts(psrhizSA.16S, function(x) x/sum(x))
filterpsrhizSA.16S.rab <- phyloseq::genefilter_sample(psrhizSA.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizSA.16Score.rab <- prune_taxa(filterpsrhizSA.16S.rab, psrhizSA.16S.rab)
rhizSA.16Score.genus.rab <- tax_glom(rhizSA.16Score.rab, taxrank = "Genus", NArm = FALSE)
#two NA are being combined. Need to make sure they stay separate
rhizSA.16Score.genus.rab <-  renameTaxa(rhizSA.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rhizSA.16Score.genus.rab.df <- psmelt(rhizSA.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSA.16Score.genus.rab.df$Genus <- as.character(rhizSA.16Score.genus.rab.df$Genus)
ggplot(rhizSA.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR  season
psrhizSR.16S.rab <- transform_sample_counts(psrhizSR.16S, function(x) x/sum(x))
filterpsrhizSR.16S.rab <- phyloseq::genefilter_sample(psrhizSR.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizSR.16Score.rab <- prune_taxa(filterpsrhizSR.16S.rab, psrhizSR.16S.rab)
rhizSR.16Score.genus.rab <- tax_glom(rhizSR.16Score.rab, taxrank = "Genus", NArm = FALSE)
#one NA are being combined. Need to make sure they stay separate
rhizSR.16Score.genus.rab <-  renameTaxa(rhizSR.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rhizSR.16Score.genus.rab.df <- psmelt(rhizSR.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSR.16Score.genus.rab.df$Genus <- as.character(rhizSR.16Score.genus.rab.df$Genus)
ggplot(rhizSR.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS  season
psrhizSS.16S.rab <- transform_sample_counts(psrhizSS.16S, function(x) x/sum(x))
filterpsrhizSS.16S.rab <- phyloseq::genefilter_sample(psrhizSS.16S.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizSS.16Score.rab <- prune_taxa(filterpsrhizSS.16S.rab, psrhizSS.16S.rab)
rhizSS.16Score.genus.rab <- tax_glom(rhizSS.16Score.rab, taxrank = "Genus", NArm = FALSE)
#one NA are being combined. Need to make sure they stay separate
rhizSS.16Score.genus.rab <-  renameTaxa(rhizSS.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rhizSS.16Score.genus.rab.df <- psmelt(rhizSS.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSS.16Score.genus.rab.df$Genus <- as.character(rhizSS.16Score.genus.rab.df$Genus)
ggplot(rhizSS.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")




# Merge Rhizosphere 16S dataframes together and make a heatmap of the 1% abundant or higher taxa ----

allrhiz.16Score.genus.df <- rbind(rhizS3.16Score.genus.rab.df, rhizSA.16Score.genus.rab.df, rhizSR.16Score.genus.rab.df, rhizSS.16Score.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allrhiz.16Score.genus.df, "allrhiz.16Score.genus.csv")
#import back in
allrhiz.16Score.genus.df <- read.csv("fix.allrhiz.16Score.genus.csv")


#WINNER
#max % is 34.08
#need to change the values depending on how you want to view it. Values are 0 to 1.  
ggplot(allrhiz.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.1, 0.2, 0.30 ,0.3408), values = c("0", "0.2", "1"))

#break


# Roots 16S taxa with 2% abundance or higher. ----
#S3  season
psrootS3.16S.rab <- transform_sample_counts(psrootS3.16S, function(x) x/sum(x))
filterpsrootS3.16S.rab <- phyloseq::genefilter_sample(psrootS3.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootS3.16Score.rab <- prune_taxa(filterpsrootS3.16S.rab, psrootS3.16S.rab)
rootS3.16Score.genus.rab <- tax_glom(rootS3.16Score.rab, taxrank = "Genus", NArm = FALSE)
#4 NA are being combined. Need to make sure they stay separate
rootS3.16Score.genus.rab <-  renameTaxa(rootS3.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rootS3.16Score.genus.rab.df <- psmelt(rootS3.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rootS3.16Score.genus.rab.df$Genus <- as.character(rootS3.16Score.genus.rab.df$Genus)
ggplot(rootS3.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SA  season
psrootSA.16S.rab <- transform_sample_counts(psrootSA.16S, function(x) x/sum(x))
filterpsrootSA.16S.rab <- phyloseq::genefilter_sample(psrootSA.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootSA.16Score.rab <- prune_taxa(filterpsrootSA.16S.rab, psrootSA.16S.rab)
rootSA.16Score.genus.rab <- tax_glom(rootSA.16Score.rab, taxrank = "Genus", NArm = FALSE)
#two NA are being combined. Need to make sure they stay separate
rootSA.16Score.genus.rab <-  renameTaxa(rootSA.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rootSA.16Score.genus.rab.df <- psmelt(rootSA.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rootSA.16Score.genus.rab.df$Genus <- as.character(rootSA.16Score.genus.rab.df$Genus)
ggplot(rootSA.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SR  season
psrootSR.16S.rab <- transform_sample_counts(psrootSR.16S, function(x) x/sum(x))
filterpsrootSR.16S.rab <- phyloseq::genefilter_sample(psrootSR.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootSR.16Score.rab <- prune_taxa(filterpsrootSR.16S.rab, psrootSR.16S.rab)
rootSR.16Score.genus.rab <- tax_glom(rootSR.16Score.rab, taxrank = "Genus", NArm = FALSE)
#4 NA are being combined. Need to make sure they stay separate
rootSR.16Score.genus.rab <-  renameTaxa(rootSR.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rootSR.16Score.genus.rab.df <- psmelt(rootSR.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rootSR.16Score.genus.rab.df$Genus <- as.character(rootSR.16Score.genus.rab.df$Genus)
ggplot(rootSR.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


#SS  season
psrootSS.16S.rab <- transform_sample_counts(psrootSS.16S, function(x) x/sum(x))
filterpsrootSS.16S.rab <- phyloseq::genefilter_sample(psrootSS.16S.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootSS.16Score.rab <- prune_taxa(filterpsrootSS.16S.rab, psrootSS.16S.rab)
rootSS.16Score.genus.rab <- tax_glom(rootSS.16Score.rab, taxrank = "Genus", NArm = FALSE)
#one NA are being combined. Need to make sure they stay separate
rootSS.16Score.genus.rab <-  renameTaxa(rootSS.16Score.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe
rootSS.16Score.genus.rab.df <- psmelt(rootSS.16Score.genus.rab)
# convert Genus to a character vector from a factor because R
rootSS.16Score.genus.rab.df$Genus <- as.character(rootSS.16Score.genus.rab.df$Genus)
ggplot(rootSS.16Score.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


# Merge root dataframes for heatmap of 1% or higher. ---- 

allroot.16Score.genus.df <- rbind(rootS3.16Score.genus.rab.df, rootSA.16Score.genus.rab.df, rootSR.16Score.genus.rab.df, rootSS.16Score.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allroot.16Score.genus.df, "allroot.16Score.genus.csv")
#import back in
allroot.16Score.genus.df <- read.csv("fix.allroot.16Score.genus.csv")


#WINNER
#max % is 73.8
#need to change the values depending on how you want to view it. Values are 0 to 1.  
ggplot(allroot.16Score.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.1, 0.2, 0.30 , 0.4, 0.5, 0.6, 0.738), values = c("0", "0.2", "1"))


#break

# PREVALENCE plots for 16S ----

#export "allcysts.16Score.genus.df" as "allcysts.16Score.genus.csv"
#in excel, make columns for # of samples present in, Prevalence (# of samples present/total samples), Median abundance
#import back in as:
cyst.all.prevalence.df <- read.csv("CystAll16SPrevalence.csv")
cyst.all.prevalence.df$Genus <- as.character(cyst.all.prevalence.df$Genus)
cyst.all.prevalence.df$Genus <- factor(cyst.all.prevalence.df$Genus, levels=unique(cyst.all.prevalence.df$Genus))

#barplot for prevalence
#order taxa in from high to low abundance
ggplot(cyst.all.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()



#export allbulk.16Score.genus.df as allbulk.ITScore.genus.csv
#in excel, make columns for # of samples present in, Prevalence (# of samples present/total samples), Median abundance
#import back in as:
bulk.all.prevalence.df <- read.csv("BulkAll16SPrevalence.csv")
bulk.all.prevalence.df$Genus <- as.character(bulk.all.prevalence.df$Genus)
bulk.all.prevalence.df$Genus <- factor(bulk.all.prevalence.df$Genus, levels=unique(bulk.all.prevalence.df$Genus))

ggplot(bulk.all.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()



#export allrhiz.16Score.genus.df as allrhiz.16Score.genus.csv
#in excel, make columns for # of samples present in, Prevalence (# of samples present/total samples), Median abundance
#import back in as:
rhiz.prevalence.df <- read.csv("RhizPrevalence.csv")
rhiz.prevalence.df$Genus <- as.character(rhiz.prevalence.df$Genus)
rhiz.prevalence.df$Genus <- factor(rhiz.prevalence.df$Genus, levels=unique(rhiz.prevalence.df$Genus))

ggplot(rhiz.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()



#export allroot.16Score.genus.df as allroot.16Score.genus.csv
#in excel, make columns for # of samples present in, Prevalence (# of samples present/total samples), Median abundance
#import back in as:
root.prevalence.df <- read.csv("RootPrevalence.csv")
root.prevalence.df$Genus <- as.character(root.prevalence.df$Genus)
root.prevalence.df$Genus <- factor(root.prevalence.df$Genus, levels=unique(root.prevalence.df$Genus))

ggplot(root.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()


#break


#Cyst Metagenome- Import Metaphlan Taxonomy files ----

#Reformat: In excel, I removed everything not labeled to the genus.

#import OTU Table
metaphlan.otutable.g <- read.csv("/workdir/eag252/metaphlan/MetaphlanBugsList_OTUTable_Genus.csv")
metaphlan.otutable.g <- data.frame(metaphlan.otutable.g, row.names = 1)
#metaphlan.otutable <- as.matrix(sapply(metaphlan.otutable, as.numeric))
metaphlan.otutable.g <- otu_table(metaphlan.otutable.g, taxa_are_rows = FALSE)

#import Tax table
metaphlan.taxtable.g <- read.csv("/workdir/eag252/metaphlan/MetaphlanBugsList_TaxaTable_Genus.csv")
metaphlan.taxtable.g <- data.frame(metaphlan.taxtable.g, row.names = 1)
metaphlan.taxtable.g <- as.matrix(metaphlan.taxtable.g)
metaphlan.taxtable.g <- tax_table(metaphlan.taxtable.g)

taxa_names(metaphlan.otutable.g)
taxa_names(metaphlan.taxtable.g)

#import metadata 
metaphlan.metadata <- read.delim("/workdir/eag252/metaphlan/MetadataMetaphlanForR.txt", header = TRUE, row.names = 1)

#Merge otu table, taxa table, metadata into phyloseq object.
psmetaphlan.g <- phyloseq(otu_table(metaphlan.otutable.g, taxa_are_rows = FALSE), sample_data(metaphlan.metadata), tax_table(metaphlan.taxtable.g))



#Split Metaphlan phyloseq object into separate Treatments ----
metaphlancystS3 <- subset_samples(psmetaphlan.g, Treatment == "S3")
metaphlancystSA <- subset_samples(psmetaphlan.g, Treatment == "SA")
metaphlancystSR <- subset_samples(psmetaphlan.g, Treatment == "SR")
metaphlancystSS <- subset_samples(psmetaphlan.g, Treatment == "SS")


# Cyst Metagenome- Metaphlan- S3 taxa .5% abundance ----
filter.metaphlancystS3 <- phyloseq::genefilter_sample(metaphlancystS3, filterfun_sample(function(x) x / sum(x) > .005))
metaphlancystS3.core <- prune_taxa(filter.metaphlancystS3, metaphlancystS3)
metaphlancystS3.core.genus <- tax_glom(metaphlancystS3.core, taxrank = "Genus", NArm = FALSE)

#make dataframe
metaphlancystS3.core.genus.df <- psmelt(metaphlancystS3.core.genus)
# convert Genus to a character vector from a factor because R
metaphlancystS3.core.genus.df$Genus <- as.character(metaphlancystS3.core.genus.df$Genus)

ggplot(metaphlancystS3.core.genus.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


# Cyst Metagenome- Metaphlan- SA taxa .5% abundance ----
filter.metaphlancystSA <- phyloseq::genefilter_sample(metaphlancystSA, filterfun_sample(function(x) x / sum(x) > .005))
metaphlancystSA.core <- prune_taxa(filter.metaphlancystSA, metaphlancystSA)
metaphlancystSA.core.genus <- tax_glom(metaphlancystSA.core, taxrank = "Genus", NArm = FALSE)

#make dataframe
metaphlancystSA.core.genus.df <- psmelt(metaphlancystSA.core.genus)
# convert Genus to a character vector from a factor because R
metaphlancystSA.core.genus.df$Genus <- as.character(metaphlancystSA.core.genus.df$Genus)

ggplot(metaphlancystSA.core.genus.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


# Cyst Metagenome- Metaphlan- SR taxa .5% abundance ----
filter.metaphlancystSR <- phyloseq::genefilter_sample(metaphlancystSR, filterfun_sample(function(x) x / sum(x) > .005))
metaphlancystSR.core <- prune_taxa(filter.metaphlancystSR, metaphlancystSR)
metaphlancystSR.core.genus <- tax_glom(metaphlancystSR.core, taxrank = "Genus", NArm = FALSE)

#make dataframe
metaphlancystSR.core.genus.df <- psmelt(metaphlancystSR.core.genus)
# convert Genus to a character vector from a factor because R
metaphlancystSR.core.genus.df$Genus <- as.character(metaphlancystSR.core.genus.df$Genus)

ggplot(metaphlancystSR.core.genus.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


# Cyst Metagenome- Metaphlan- SS taxa .5% abundance ----
filter.metaphlancystSS <- phyloseq::genefilter_sample(metaphlancystSS, filterfun_sample(function(x) x / sum(x) > .005))
metaphlancystSS.core <- prune_taxa(filter.metaphlancystSS, metaphlancystSS)
metaphlancystSS.core.genus <- tax_glom(metaphlancystSS.core, taxrank = "Genus", NArm = FALSE)

#make dataframe
metaphlancystSS.core.genus.df <- psmelt(metaphlancystSS.core.genus)
# convert Genus to a character vector from a factor because R
metaphlancystSS.core.genus.df$Genus <- as.character(metaphlancystSS.core.genus.df$Genus)

ggplot(metaphlancystSS.core.genus.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")



#break


# Cyst Metagenome- Metaphlan Heatmap ---- 

#merge each treatments dataframe together. 
metaphlancysts.allcore.df <- rbind(metaphlancystS3.core.genus.df, metaphlancystSA.core.genus.df, metaphlancystSR.core.genus.df, metaphlancystSS.core.genus.df)

#Reorder to make mid season show up first
metaphlancysts.allcore.df$Season <- factor(metaphlancysts.allcore.df$Season, levels = c("Mid", "Fall"))

ggplot(metaphlancysts.allcore.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Treatment*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue"), breaks = c(0.0, 0.2, 0.4, 0.568), values = c("0", "0.3", "1"))

#break

# NetCoMi- cysts ----
#both seasons
filterpscyst.16S.rab <- phyloseq::genefilter_sample(pscyst.rab, filterfun_sample(function(x) x / sum(x) > 0.02))
pscyst.16S.rab.filt <- prune_taxa(filterpscyst.16S.rab, pscyst.rab)

pscyst.16S.rab.glom <- tax_glom(pscyst.16S.rab.filt, taxrank = "Genus", NArm = FALSE)

pscyst.16S.rab.rename <-  renameTaxa(pscyst.16S.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_pscystglom_spiec <- netConstruct(pscyst.16S.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_pscystglom_spiec <- netAnalyze(net_pscystglom_spiec)
summary(analyze_pscystglom_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_pscystglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Cysts glom >2% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as "NetCoMi Cysts 16S Both Seasons 2percent.pdf"


#break

# NetCoMi- bulk soil ----
# both seasons
filterpsbulk.16S.rab <- phyloseq::genefilter_sample(psbulk.rab, filterfun_sample(function(x) x / sum(x) > 0.01))
psbulk.16S.rab.filt <- prune_taxa(filterpsbulk.16S.rab, psbulk.rab)

psbulk.16S.rab.glom <- tax_glom(psbulk.16S.rab.filt, taxrank = "Genus", NArm = FALSE)

psbulk.16S.rab.rename <-  renameTaxa(psbulk.16S.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psbulkglom_spiec <- netConstruct(psbulk.16S.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_psbulkglom_spiec <- netAnalyze(net_psbulkglom_spiec)
summary(analyze_psbulkglom_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_psbulkglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Bulk glom >1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as NetCoMi Bulk 16S Both Seasons 1percent.pdf

#break

# NetCoMi- rhizosphere- No network at 1% abundance ----
#No network at 1% abundance
filterpsrhiz.16S.rab <- phyloseq::genefilter_sample(psrhiz.rab, filterfun_sample(function(x) x / sum(x) > 0.01))
psrhiz.16S.rab.filt <- prune_taxa(filterpsrhiz.16S.rab, psrhiz.rab)

psrhiz.16S.rab.glom <- tax_glom(psrhiz.16S.rab.filt, taxrank = "Genus", NArm = FALSE)

psrhiz.16S.rab.rename <-  renameTaxa(psrhiz.16S.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psrhizglom_spiec <- netConstruct(psrhiz.16S.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_psrhizglom_spiec <- netAnalyze(net_psrhizglom_spiec)
summary(analyze_psrhizglom_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_psrhizglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Rhiz glom >1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#break


# NetCoMi- roots ----
filterpsroot.16S.rab <- phyloseq::genefilter_sample(psroot.rab, filterfun_sample(function(x) x / sum(x) > 0.02))
psroot.16S.rab.filt <- prune_taxa(filterpsroot.16S.rab, psroot.rab)

psroot.16S.rab.glom <- tax_glom(psroot.16S.rab.filt, taxrank = "Genus", NArm = FALSE)

psroot.16S.rab.rename <-  renameTaxa(psroot.16S.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psrootglom_spiec <- netConstruct(psroot.16S.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_psrootglom_spiec <- netAnalyze(net_psrootglom_spiec)
summary(analyze_psrootglom_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_psrootglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Roots glom >2% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as NetCoMi Roots 16S Both Seasons 2percent.pdf

#break

#NetCoMi- by treatment- all samples have the same cut off of 1%----
ps4.rab <- transform_sample_counts(ps4, function(x) x/sum(x))

filterps4.rab <- phyloseq::genefilter_sample(ps4.rab, filterfun_sample(function(x) x / sum(x) > .01))
ps4.filt.rab <- prune_taxa(filterps4.rab, ps4.rab)

#split into treatments 

ps4.filt.S3 <- subset_samples(ps4.filt.rab, Plot == "S3")

ps4.filt.SA <- subset_samples(ps4.filt.rab, Plot == "Sa")

ps4.filt.SR <- subset_samples(ps4.filt.rab, Plot == "Sr")

ps4.filt.SS <- subset_samples(ps4.filt.rab, Plot == "Ss")

#break 

#NetCoMi- by treatment-all samples 1% abund S3----
ps4.filt.S3.glom <- tax_glom(ps4.filt.S3, taxrank = "Genus", NArm = FALSE)

ps4.filt.S3.rename <-  renameTaxa(ps4.filt.S3.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps4.S3_spiec <- netConstruct(ps4.filt.S3.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_ps4.S3_spiec <- netAnalyze(net_ps4.S3_spiec)
summary(analyze_ps4.S3_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_ps4.S3_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "S3 16S all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as "NetCoMi S3 16S all samples 1percent abund.pdf"

#break

#NetCoMi- by treatment-all samples 1% abund SA----
ps4.filt.SA.glom <- tax_glom(ps4.filt.SA, taxrank = "Genus", NArm = FALSE)

ps4.filt.SA.rename <-  renameTaxa(ps4.filt.SA.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps4.SA_spiec <- netConstruct(ps4.filt.SA.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_ps4.SA_spiec <- netAnalyze(net_ps4.SA_spiec)
summary(analyze_ps4.SA_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_ps4.SA_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "SA 16S all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as "NetCoMi SA 16S all samples 1percent abund.pdf"

#break

#NetCoMi- by treatment-all samples 1% abund SR----
ps4.filt.SR.glom <- tax_glom(ps4.filt.SR, taxrank = "Genus", NArm = FALSE)

ps4.filt.SR.rename <-  renameTaxa(ps4.filt.SR.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps4.SR_spiec <- netConstruct(ps4.filt.SR.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_ps4.SR_spiec <- netAnalyze(net_ps4.SR_spiec)
summary(analyze_ps4.SR_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_ps4.SR_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "SR 16S all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as "NetCoMi SR 16S all samples 1percent abund.pdf"

#break

#NetCoMi- by treatment-all samples 1% abund SS----
ps4.filt.SS.glom <- tax_glom(ps4.filt.SS, taxrank = "Genus", NArm = FALSE)

ps4.filt.SS.rename <-  renameTaxa(ps4.filt.SS.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps4.SS_spiec <- netConstruct(ps4.filt.SS.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_ps4.SS_spiec <- netAnalyze(net_ps4.SS_spiec)
summary(analyze_ps4.SS_spiec)
#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 
plot(analyze_ps4.SS_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "SS 16S all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)
#Saved as "NetCoMi SS 16S all samples 1percent abund.pdf"

#break

# ANCOM-BC- making a dataset with the same thresholds as the top taxa ----

#Make a phyloseq object with the same taxa that are included in abundance thresholds. 
#Cysts = 2% or higher, Soils = 1% or higher, Roots = 2% or higher. 

pscyst.rab <- transform_sample_counts(pscyst, function(x) x/sum(x))
filter.pscyst.rab <- phyloseq::genefilter_sample(pscyst.rab, filterfun_sample(function(x) x / sum(x) > .02))
#pull from the non relative abundance bulk file
forancom.filter.pscyst <- prune_taxa(filter.pscyst.rab, pscyst)

pscystmid.rab <- transform_sample_counts(pscystmid, function(x) x/sum(x))
filter.pscystmid.rab <- phyloseq::genefilter_sample(pscystmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.pscystmid <- prune_taxa(filter.pscystmid.rab, pscystmid)


pscystfall.rab <- transform_sample_counts(pscystfall, function(x) x/sum(x))
filter.pscystfall.rab <- phyloseq::genefilter_sample(pscystfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.pscystfall <- prune_taxa(filter.pscystfall.rab, pscystfall)


psbulk.rab <- transform_sample_counts(psbulk, function(x) x/sum(x))
filter.psbulk.rab <- phyloseq::genefilter_sample(psbulk.rab, filterfun_sample(function(x) x / sum(x) > .01))
#pull from the non relative abundance bulk file
forancom.filter.psbulk <- prune_taxa(filter.psbulk.rab, psbulk)


psbulkmid.rab <- transform_sample_counts(psbulkmid, function(x) x/sum(x))
filter.psbulkmid.rab <- phyloseq::genefilter_sample(psbulkmid.rab, filterfun_sample(function(x) x / sum(x) > .01))
forancom.filter.psbulkmid <- prune_taxa(filter.psbulkmid.rab, psbulkmid)


psbulkfall.rab <- transform_sample_counts(psbulkfall, function(x) x/sum(x))
filter.psbulkfall.rab <- phyloseq::genefilter_sample(psbulkfall.rab, filterfun_sample(function(x) x / sum(x) > .01))
forancom.filter.psbulkfall <- prune_taxa(filter.psbulkfall.rab, psbulkfall)


psrhiz.rab <- transform_sample_counts(psrhiz, function(x) x/sum(x))
filter.psrhiz.rab <- phyloseq::genefilter_sample(psrhiz.rab, filterfun_sample(function(x) x / sum(x) > .01))
forancom.filter.psrhiz <- prune_taxa(filter.psrhiz.rab, psrhiz)

psroot.rab <- transform_sample_counts(psroot, function(x) x/sum(x))
filter.psroot.rab <- phyloseq::genefilter_sample(psroot.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.psroot <- prune_taxa(filter.psroot.rab, psroot)

#Ps5 has abundance cut offs for each sample! 
ps5 <- merge_phyloseq(forancom.filter.pscyst, forancom.filter.psbulk, forancom.filter.psrhiz, forancom.filter.psroot)
#check to see if it works
testancombc.ps5 <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                            group = "Sample_Type", n_cl = 8, global = TRUE)

#Bulk as comparison
sample_data(ps5)$Sample_Type <- as.factor(sample_data(ps5)$Sample_Type)
ancombc.bulk <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.bulk = ancombc.bulk$res
res.bulk_global = ancombc.bulk$res_global


#Cyst 
sample_data(ps5)$Sample_Type <- relevel(sample_data(ps5)$Sample_Type, "Cyst")
ancombc.cyst <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.cyst = ancombc.cyst$res
res.cyst_global = ancombc.cyst$res_global


#Rhizosphere
sample_data(ps5)$Sample_Type <- relevel(sample_data(ps5)$Sample_Type, "Rhizosphere")
ancombc.rhiz <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.rhiz = ancombc.rhiz$res
res.rhiz_global = ancombc.rhiz$res_global


#Root
sample_data(ps5)$Sample_Type <- relevel(sample_data(ps5)$Sample_Type, "Roots")
ancombc.root <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.root = ancombc.root$res
res.root_global = ancombc.root$res_global

write.csv(res.bulk, "res.bulk.csv")
write.csv(res.cyst, "res.cyst.csv")
write.csv(res.rhiz, "res.rhiz.csv")
write.csv(res.root, "res.root.csv")

#Using the res.cyst.csv 
#Removed anything that was false across all treatments and saved as filt.rest.cyst.csv
#made 3 csv files that have the 4 columns. Genera, lfc T or F values, lfc number values, and the Cyst vs a sample type
#I can merge these files now to make a dataframe and try a heatmap? 

cvb <- read.csv("cystvsbulk.csv")
cvrh <- read.csv("cystvsrhiz.csv")
cvr <- read.csv("cystvsroot.csv")

cystvsall.ancom.df <- rbind(cvb, cvrh, cvr)

ggplot(cystvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Cysts vs All Ancom.pdf"

#bulk
bvc <- read.csv("bulkvscyst.csv")
bvrh <- read.csv("bulkvsrhiz.csv")
bvr <- read.csv("bulkvsroot.csv")

bulkvsall.ancom.df <- rbind(bvc, bvrh, bvr)

ggplot(bulkvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Bulk vs All Ancom.pdf"

#rhiz
rhvb <- read.csv("rhizvsbulk.csv")
rhvc <- read.csv("rhizvscyst.csv")
rhvr <- read.csv("rhizvsroot.csv")

rhizvsall.ancom.df <- rbind(rhvb, rhvc, rhvr)

ggplot(rhizvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Rhiz vs All Ancom.pdf"

#roots
rvb <- read.csv("rootvsbulk.csv")
rvc <- read.csv("rootvscyst.csv")
rvrh <- read.csv("rootvsrhiz.csv")

rootvsall.ancom.df <- rbind(rvb, rvc, rvrh)

ggplot(rootvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Roots vs All Ancom.pdf"

#all sample types & comparisons in one.
all.ancom.df <- rbind(cystvsall.ancom.df, bulkvsall.ancom.df, rhizvsall.ancom.df, rootvsall.ancom.df)

ggplot(all.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "All Ancom in one plot.pdf"

#break

#ANCOM-BC- Entire dataset, compare treatments- nothing was significant ---- 

#S3
sample_data(ps5)$Plot <- as.factor(sample_data(ps5)$Plot)
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "S3")
ancombc.ps5.S3 <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.ps5.S3 = ancombc.ps5.S3$res


#Sa
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "Sa")
ancombc.ps5.Sa <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.ps5.Sa = ancombc.ps5.Sa$res


#Sr
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "Sr")
ancombc.ps5.Sr <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.ps5.Sr = ancombc.ps5.Sr$res


#Ss
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "Ss")
ancombc.ps5.Ss <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.ps5.Ss = ancombc.ps5.Ss$res

write.csv(res.ps5.S3, "res.ps5.S3.csv")
write.csv(res.ps5.Sa, "res.ps5.SA.csv")
write.csv(res.ps5.Sr, "res.ps5.SR.csv")
write.csv(res.ps5.Ss, "res.ps5.SS.csv")

#break


#ANCOM-BC- Just cysts- compare the treatments ---- 
#Just cysts- compare treatments

#using the filtered phyloseq object of cysts
forancom.filter.pscyst

sample_data(forancom.filter.pscyst)$Plot <- as.factor(sample_data(forancom.filter.pscyst)$Plot)
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "S3")
ancombc.cystS3 <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.cystS3 = ancombc.cystS3$res
write.csv(res.cystS3, "res.cystS3.csv")

#Sa
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "Sa")
ancombc.cystSA <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.cystSA = ancombc.cystSA$res
write.csv(res.cystSA, "res.cystSA.csv")

#Sr
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "Sr")
ancombc.cystSR <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.cystSR = ancombc.cystSR$res
write.csv(res.cystSR, "res.cystSR.csv")

#Ss
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "Ss")
ancombc.cystSS <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.cystSS = ancombc.cystSS$res
write.csv(res.cystSS, "res.cystSS.csv")


#Removed anything that was false across all treatments and saved as filt.rest.cyst.csv
#made new csv files that have the 4 columns: Genera, lfc T or F values, lfc number values, and the Cyst vs a sample type
#I can merge these files now to make a dataframe and try a heatmap 

cyst.savs3 <- read.csv("cystSAvcystS3.csv")
cyst.savsr <- read.csv("cystSAvcystSR.csv")
cyst.savss <- read.csv("cystSAvcystSS.csv")

cyst.srvs3 <- read.csv("cystSRvcystS3.csv")
cyst.srvsa <- read.csv("cystSRvcystSA.csv")
cyst.srvss <- read.csv("cystSRvcystSS.csv")

cyst.ssvs3 <- read.csv("cystSSvcystS3.csv")
cyst.ssvsa <- read.csv("cystSSvcystSA.csv")
cyst.ssvsr <- read.csv("cystSSvcystSR.csv")

ancom.cyst.savall.df <- rbind(cyst.savs3, cyst.savsr, cyst.savss)

ancom.cyst.srvall.df <- rbind(cyst.srvs3, cyst.srvsa, cyst.srvss)

ancom.cyst.ssvall.df <- rbind(cyst.ssvs3, cyst.ssvsa, cyst.ssvsr)


ancom.cyst.df <- rbind(ancom.cyst.savall.df, ancom.cyst.srvall.df, ancom.cyst.ssvall.df)

ggplot(ancom.cyst.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom Cyst- merged treatments.pdf"

#break


#ANCOM-BC- Just mid season cysts- compare the treatments ---- 
#using the filtered phyloseq object of cysts
forancom.filter.pscystmid 

sample_data(forancom.filter.pscystmid)$Plot <- as.factor(sample_data(forancom.filter.pscystmid)$Plot)
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "S3")
ancombc.cystmidS3 <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.cystmidS3 = ancombc.cystmidS3$res
#write.csv(res.cystmidS3, "res.cystmidS3.csv")
#Copied from dataframe into an excel sheet to save time.

#Sa
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "Sa")
ancombc.cystmidSA <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)
res.cystmidSA = ancombc.cystmidSA$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystmidSA, "res.cystmidSA.csv")

#Sr
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "Sr")
ancombc.cystmidSR <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)
res.cystmidSR = ancombc.cystmidSR$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystmidSR, "res.cystmidSR.csv")

#Ss
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "Ss")
ancombc.cystmidSS <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.cystmidSS = ancombc.cystmidSS$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystmidSS, "res.cystmidSS.csv")

#Removed anything that was false across all treatments
#made new csv files that have the 4 columns: Genera, lfc T or F values, lfc number values, and the sample vs sample 
#I can merge these files now to make a dataframe and try a heatmap 

cystmid.s3vsa <- read.csv("cystmidS3vcystmidSA.csv")
cystmid.s3vsr <- read.csv("cystmidS3vcystmidSR.csv")

cystmid.savs3 <- read.csv("cystmidSAvcystmidS3.csv")
cystmid.savsr <- read.csv("cystmidSAvcystmidSR.csv")
cystmid.savss <- read.csv("cystmidSAvcystmidSS.csv")

cystmid.srvs3 <- read.csv("cystmidSRvcystmidS3.csv")
cystmid.srvsa <- read.csv("cystmidSRvcystmidSA.csv")

cystmid.ssvsa <- read.csv("cystmidSSvcystmidSA.csv")


ancom.cystmid.s3vall.df <- rbind(cystmid.s3vsa, cystmid.s3vsr)

ancom.cystmid.savall.df <- rbind(cystmid.savs3, cystmid.savsr, cystmid.savss)

ancom.cystmid.srvall.df <- rbind(cystmid.srvs3, cystmid.srvsa)


ancom.cystmid.df <- rbind(ancom.cystmid.s3vall.df, ancom.cystmid.savall.df, ancom.cystmid.srvall.df, cystmid.ssvsa)

ggplot(ancom.cystmid.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom Cyst Mid- merge treatments.pdf"

#break


#ANCOM-BC- Just fall season cysts- compare the treatments ---- 
#using the filtered phyloseq object of cysts
forancom.filter.pscystfall 

sample_data(forancom.filter.pscystfall)$Plot <- as.factor(sample_data(forancom.filter.pscystfall)$Plot)
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "S3")
ancombc.cystfallS3 <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.cystfallS3 = ancombc.cystfallS3$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystfallS3, "res.cystfallS3.csv")

#Sa
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "Sa")
ancombc.cystfallSA <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)
res.cystfallSA = ancombc.cystfallSA$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystfallSA, "res.cystfallSA.csv")

#Sr
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "Sr")
ancombc.cystfallSR <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)
res.cystfallSR = ancombc.cystfallSR$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystfallSR, "res.cystfallSR.csv")

#Ss
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "Ss")
ancombc.cystfallSS <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.cystfallSS = ancombc.cystfallSS$res
#Copied from dataframe into an excel sheet to save time.  write.csv(res.cystfallSS, "res.cystfallSS.csv")


#Removed anything that was false across all treatments
#made new csv files that have the 4 columns: Genera, lfc T or F values, lfc number values, and the sample vs sample 
#I can merge these files now to make a dataframe and try a heatmap 

cystfall.savsr <- read.csv("cystfallSAvcystfallSR.csv")
cystfall.savss <- read.csv("cystfallSAvcystfallSS.csv")

cystfall.srvsa <- read.csv("cystfallSRvcystfallSA.csv")

cystfall.ssvsa <- read.csv("cystfallSSvcystfallSA.csv")


ancom.cystfall.df <- rbind(cystfallSAvcystfallSR.csv, cystfallSAvcystfallSS.csv, ancom.cystfall.savall.df, cystfall.srvsa, cystfall.ssvsa)

ggplot(ancom.cystfall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom Cyst Fall- merged treatments.pdf"

#break

#ANCOM-BC- Just Cyst- compare seasons ---- 

ancombc.cyst.season <- ancombc2(forancom.filter.pscyst, assay_name = "counts", tax_level = "Genus", fix_formula = "Season",
                                group = "Season", n_cl = 8, global = TRUE)

res.cyst.season = ancombc.cyst.season$res
write.csv(res.cyst.season, "res.cyst.season.csv")

#Removed anything that was false across all treatments
#made new csv files that have the 4 columns: Genera, lfc T or F values, lfc number values, and the sample vs sample 
#I can merge these files now to make a dataframe and try a heatmap 

ancom.cyst.season <- read.csv("cyst.season.csv")

ggplot(ancom.cyst.season, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Cysts- seasons.pdf"

#break

#ANCOM-BC- Just bulk- compare the treatments- Nothing was significant! ---- 
forancom.filter.psbulk

sample_data(forancom.filter.psbulk)$Plot <- as.factor(sample_data(forancom.filter.psbulk)$Plot)
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "S3")
ancombc.bulkS3 <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.bulkS3 = ancombc.bulkS3$res
#write.csv(res.bulkS3, "res.bulkS3.csv")

#Sa
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "Sa")
ancombc.bulkSA <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.bulkSA = ancombc.bulkSA$res
#write.csv(res.bulkSA, "res.bulkSA.csv")

#Sr
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "Sr")
ancombc.bulkSR <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.bulkSR = ancombc.bulkSR$res
#write.csv(res.bulkSR, "res.bulkSR.csv")

#Ss
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "Ss")
ancombc.bulkSS <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.bulkSS = ancombc.bulkSS$res
#write.csv(res.bulkSS, "res.bulkSS.csv")

#Nothing was significant!

#break

#ANCOM-BC- Just mid bulk- compare the treatments ----
forancom.filter.psbulkmid <- subset_samples(forancom.filter.psbulk, Season == "Mid")

sample_data(forancom.filter.psbulkmid)$Plot <- as.factor(sample_data(forancom.filter.psbulkmid)$Plot)
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "S3")
ancombc.bulkmidS3 <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.bulkmidS3 = ancombc.bulkmidS3$res
#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidS3, "res.bulkmidS3.csv")

#Sa
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "Sa")
ancombc.bulkmidSA <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)
res.bulkmidSA = ancombc.bulkmidSA$res
#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidSA, "res.bulkmidSA.csv")

#Sr
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "Sr")
ancombc.bulkmidSR <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)
res.bulkmidSR = ancombc.bulkmidSR$res
#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidSR, "res.bulkmidSR.csv")

#Ss
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "Ss")
ancombc.bulkmidSS <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.bulkmidSS = ancombc.bulkmidSS$res
#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidSS, "res.bulkmidSS.csv")


bulkmid.s3vsa <- read.csv("bulkmidS3vbulkmidSA.csv")

bulkmid.savs3 <- read.csv("bulkmidSAvbulkmidS3.csv")
bulkmid.savsr <- read.csv("bulkmidSAvbulkmidSR.csv")
bulkmid.savss <- read.csv("bulkmidSAvbulkmidSS.csv")

bulkmid.srvsa <- read.csv("bulkmidSRvbulkmidSA.csv")

bulkmid.ssvsa <- read.csv("bulkmidSSvbulkmidSA.csv")


ggplot(bulkmid.s3vsa, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Bulk Mid S3 vs all.pdf"

ancom.bulkmid.savall.df <- rbind(bulkmid.savs3, bulkmid.savsr, bulkmid.savss)

ggplot(ancom.bulkmid.savall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Bulk Mid SA vs all.pdf"

ggplot(bulkmid.srvsa, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Bulk Mid SR vs all.pdf"

ggplot(bulkmid.ssvsa, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Bulk Mid SS vs all.pdf"


ancom.bulkmid.df <- rbind(bulkmid.s3vsa, ancom.bulkmid.savall.df, bulkmid.srvsa, bulkmid.ssvsa)

ggplot(ancom.bulkmid.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom Bulk Mid- merged treatments.pdf"

#break 


#ANCOM-BC- Just fall bulk- compare the treatments NONE ARE SIGNIFICANT ----
forancom.filter.psbulkfall <- subset_samples(forancom.filter.psbulk, Season == "Fall")

sample_data(forancom.filter.psbulkfall)$Plot <- as.factor(sample_data(forancom.filter.psbulkfall)$Plot)
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "S3")
ancombc.bulkfallS3 <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.bulkfallS3 = ancombc.bulkfallS3$res
#write.csv(res.bulkfallS3, "res.bulkfallS3.csv")

#Sa
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "Sa")
ancombc.bulkfallSA <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)
res.bulkfallSA = ancombc.bulkfallSA$res
#write.csv(res.bulkfallSA, "res.bulkfallSA.csv")

#Sr
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "Sr")
ancombc.bulkfallSR <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)
res.bulkfallSR = ancombc.bulkfallSR$res
#write.csv(res.bulkfallSR, "res.bulkfallSR.csv")

#Ss
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "Ss")
ancombc.bulkfallSS <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.bulkfallSS = ancombc.bulkfallSS$res
#write.csv(res.bulkfallSS, "res.bulkfallSS.csv")

#break

#ANCOM-BC- Just Bulk- compare seasons ---- 

ancombc.bulk.season <- ancombc2(forancom.filter.psbulk, assay_name = "counts", tax_level = "Genus", fix_formula = "Season",
                                group = "Season", n_cl = 8, global = TRUE)

res.bulk.season = ancombc.bulk.season$res
write.csv(res.bulk.season, "res.bulk.season.csv")

ancom.bulk.season <- read.csv("bulk.season.csv")

ggplot(ancom.bulk.season, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Bulk- seasons.pdf"

#break

#ANCOM-BC- Just rhiz- compare the treatments ---- 
forancom.filter.psrhiz

sample_data(forancom.filter.psrhiz)$Plot <- as.factor(sample_data(forancom.filter.psrhiz)$Plot)
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "S3")
ancombc.rhizS3 <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rhizS3 = ancombc.rhizS3$res
write.csv(res.rhizS3, "res.rhizS3.csv")

#Sa
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "Sa")
ancombc.rhizSA <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.rhizSA = ancombc.rhizSA$res
#write.csv(res.rhizSA, "res.rhizSA.csv")

#Sr
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "Sr")
ancombc.rhizSR <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.rhizSR = ancombc.rhizSR$res
write.csv(res.rhizSR, "res.rhizSR.csv")

#Ss
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "Ss")
ancombc.rhizSS <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rhizSS = ancombc.rhizSS$res
write.csv(res.rhizSS, "res.rhizSS.csv")

rhiz.s3vsa <- read.csv("rhizS3vrhizSA.csv")
rhiz.s3vsr <- read.csv("rhizS3vrhizSR.csv")
rhiz.s3vss <- read.csv("rhizS3vrhizSS.csv")

rhiz.srvs3 <- read.csv("rhizSRvrhizS3.csv")
rhiz.srvsa <- read.csv("rhizSRvrhizSA.csv")
rhiz.srvss <- read.csv("rhizSRvrhizSS.csv")

rhiz.ssvs3 <- read.csv("rhizSSvrhizS3.csv")
rhiz.ssvsa <- read.csv("rhizSSvrhizSA.csv")
rhiz.ssvsr <- read.csv("rhizSSvrhizSR.csv")


ancom.rhiz.s3vall.df <- rbind(rhiz.s3vsa, rhiz.s3vsr, rhiz.s3vss)

ggplot(ancom.rhiz.s3vall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom-Rhiz S3 vs All.pdf"

ancom.rhiz.srvall.df <- rbind(rhiz.srvs3, rhiz.srvsa, rhiz.srvss)

ggplot(ancom.rhiz.srvall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Rhiz SR vs All.pdf"

ancom.rhiz.ssvall.df <- rbind(rhiz.ssvs3, rhiz.ssvsa, rhiz.ssvsr)

ggplot(ancom.rhiz.ssvall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom- Rhiz SS vs All.pdf"

ancom.rhiz.df <- rbind(ancom.rhiz.s3vall.df, ancom.rhiz.srvall.df, ancom.rhiz.ssvall.df)

ggplot(ancom.rhiz.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom Rhiz- merged treatments.pdf"

#break

#ANCOM-BC- Just roots- compare the treatments-- nothing was significant---- 
forancom.filter.psroot

sample_data(forancom.filter.psroot)$Plot <- as.factor(sample_data(forancom.filter.psroot)$Plot)
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "S3")
ancombc.rootS3 <- ancombc2(forancom.filter.psroot ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rootS3 = ancombc.rootS3$res
write.csv(res.rootS3, "res.rootsS3.csv")

#Sa
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "Sa")
ancombc.rootSA <- ancombc2(forancom.filter.psroot ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.rootSA = ancombc.rootSA$res
write.csv(res.rootSA, "res.rootsSA.csv")

#Sr
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "Sr")
ancombc.rootSR <- ancombc2(forancom.filter.psroot, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.rootSR = ancombc.rootSR$res
write.csv(res.rootSR, "res.rootSR.csv")

#Ss
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "Ss")
ancombc.rootSS <- ancombc2(forancom.filter.psroot ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rootSS = ancombc.rootSS$res
write.csv(res.rootsSS, "res.rootsSS.csv")

#nothing was significant! 
#break
