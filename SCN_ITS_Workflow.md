# Emily A Green PhD. 

### Cyst samples that are not present in the cyst metagenome were removed
### All S1 samples removed
### This code reflects being run in the /workdir/eag252/LTR_ITS_S3ARS/ directory.
### If being run in the "storage" directory, file paths will likely need to be updated. 

## dada2 processing
library(dada2)

path <- "/workdir/eag252/LTR_ITS_S3ARS/"

list.files(path)

#### Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
#### Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[20:40])

#The negative controls are pretty bad that is fine!!

#There are many more reads than previous sequencing. The last 100 bp or so don't look good though

plotQualityProfile(fnRs[20:40])

#Reverse reads are pretty bad

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names

names(filtRs) <- sample.names

#ITS pipeline. Removes reads less than 50bp

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

#Some input samples had no reads pass the filter.

head(out)

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#### Inspect the merger data.frame from the first sample

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

#### Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

##### If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

#Export file to see all read counts in dada2 pipeline

write.csv(track, "allITScounts.csv")


### UNITE fungal database from https://doi.plutof.ut.ee/doi/10.15156/BIO/2483911
#### sh_general_release_27.10.2022.tgz
taxa <- assignTaxonomy(seqtab.nochim, "/workdir/eag252/LTR_ITS_S3ARS/sh_general_release_dynamic_27.10.2022.fasta", multithread=TRUE)

taxa.print <- taxa 

rownames(taxa.print) <- NULL

head(taxa.print)


## Load packages
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


## Read in metadata
meta <- read.csv("/workdir/eag252/LTR_ITS_S3ARS/LTRITS_S3ARSMetadata.csv",header = TRUE, row.names = 1)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(meta), tax_table(taxa))

sample_data(ps)

ps

## Filter Taxa
### Remove any taxa not identified to the phylum level. 
sort(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))

sort(get_taxa_unique(ps, taxonomic.rank = "Phylum"))

View(tax_table(ps))

ps1 <- subset_taxa(ps, !Phylum == "NA")

View(tax_table(ps1))

## Decontam
#### Checks for potential contaminates that are in the negative controls & present in the samples
library(decontam)

sample_variables(ps1)

df <- as.data.frame(sample_data(ps1))

df$LibrarySize <- sample_sums(ps1)

df <- df[order(df$LibrarySize),]

df$Index <- seq(nrow(df))

ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(ps1)$is.neg <- sample_data(ps1)$Sample_or_Control == "Control Sample"

contamdf.prev <- isContaminant(ps1, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant)

sample_sums(ps1)

View(sample_sums(ps1))

## Quick check on alpha diversity 
plot_richness(ps1, x="SampleType", measures=c("Shannon", "Simpson"), color="Sample_or_Control")

## Remove samples with less than 100 reads
ps2 <- prune_samples(sample_sums(ps1) > 100, ps1)
View(sample_sums(ps2))
View(sample_data(ps2))

## Beta Diversity- Unifrac- all samples, negative controls & mock community 
ps2t <- ps2
#Unifrac
ps2t.tree <- rtree(ntaxa(ps2t), rooted = TRUE, tip.label = taxa_names(ps2t))
#plot(ps2t.tree)
ps2t <- merge_phyloseq(ps2t, ps2t.tree)

#Unweighted
ps2tunifrac_dist <- phyloseq::distance(ps2t, method="unifrac", weighted=F)
ordination = ordinate(ps2t, method="PCoA", distance=ps2tunifrac_dist)
plot_ordination(ps2t, ordination, color = "Sample_or_Control", shape = "Sample_Type") + theme(aspect.ratio=1)

#statistics
#permanova
ps2t.tree.data <- data.frame(sample_data(ps2t))
adonis2(ps2tunifrac_dist ~SampleType*Sample_or_Control, ps2t.tree.data)

#Weighted
ps2tunifrac_dist <- phyloseq::distance(ps2t, method="unifrac", weighted=T)
ordination = ordinate(ps2t, method="PCoA", distance=ps2tunifrac_dist)
plot_ordination(ps2t, ordination, color = "Sample_or_Control", shape = "Sample_Type") + theme(aspect.ratio=1)

#statistics
#permanova
ps2t.tree.data <- data.frame(sample_data(ps2t))
adonis2(ps2tunifrac_dist ~SampleType*Sample_or_Control, ps2t.tree.data)



#Remove negative controls & mock communities 
ps3 <- subset_samples(ps2, Sample_or_Control == "Sample")
View(sample_data(ps3))

#Beta Diversity- Unifrac- All samples with no negative controls/mock community 
ps3t <- ps3
#Unifrac
ps3t.tree <- rtree(ntaxa(ps3t), rooted = TRUE, tip.label = taxa_names(ps3t))
#plot(ps3t.tree)
ps3t <- merge_phyloseq(ps3t, ps3t.tree)

#Unweighted
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=F)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "Sample_Type", shape = "Season") + theme(aspect.ratio=1
)+ geom_point(size = 2)+ scale_color_manual(values = c("Bulk soil" = "#F8766D", "Roots" = "#7CAE00", "Cyst" = "#C77CFF", "Rhizosphere" = "#00BFC4"))

#statistics
#permanova
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~Sample_Type*Season*Plot, ps3t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(ps3tunifrac_dist, ps3t.tree.data$Sample_Type, p.adjust.m = "holm", perm = 9999)

pairwiseAdonis::pairwise.adonis(ps3tunifrac_dist, ps3t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
ps3tu.bdisper <- betadisper(ps3tunifrac_dist, ps3t.tree.data$Sample_Type, "centroid", bias.adjust = TRUE)
anova(ps3tu.bdisper)
permutest(ps3tu.bdisper, pairwise = TRUE)



#Weighted
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=T)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "Sample_Type", shape = "Season") + theme(aspect.ratio=1
)+ geom_point(size = 2) + scale_color_manual(values = c("Bulk soil" = "#F8766D", "Roots" = "#7CAE00", "Cyst" = "#C77CFF", "Rhizosphere" = "#00BFC4"))

#statistics
#permanova
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~Sample_Type*Season*Plot, ps3t.tree.data)

#permanova posthoc
pairwiseAdonis::pairwise.adonis(ps3tunifrac_dist, ps3t.tree.data$Sample_Type, p.adjust.m = "holm", perm = 9999)

pairwiseAdonis::pairwise.adonis(ps3tunifrac_dist, ps3t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
ps3tu.bdisper <- betadisper(ps3tunifrac_dist, ps3t.tree.data$Sample_Type, "centroid", bias.adjust = TRUE)
anova(ps3tu.bdisper)
permutest(ps3tu.bdisper, pairwise = TRUE)


#Split all samples into separate phyloseq objects
pssoil <- subset_samples(ps3, SampleType == "Soil")
psbulk <- subset_samples(ps3, SoilType == "Bulk soil")
psrhiz <- subset_samples(ps3, SoilType == "Rhizosphere")
pscyst <- subset_samples(ps3, SampleType == "Cyst")
psroot <- subset_samples(ps3, SampleType == "Roots")

View(sample_data(pssoil))
View(tax_table(pssoil))

View(sample_data(pscyst))
View(tax_table(pscyst))

View(sample_data(psroot))
View(tax_table(psroot))


#Remove taxa with 0 reads from each of the individual sample type phyloseq objects----
#4500 taxa
pssoil <- prune_taxa(taxa_sums(pssoil) > 0, pssoil)
View(taxa_sums(pssoil))
View(tax_table(pssoil))
#removed ~500 taxa

psbulk <- prune_taxa(taxa_sums(psbulk) > 0, psbulk)
View(tax_table(psbulk))
#remove 1055

psrhiz <- prune_taxa(taxa_sums(psrhiz) > 0, psrhiz)
View(tax_table(psrhiz))
#removed 2821

pscyst <- prune_taxa(taxa_sums(pscyst) > 0, pscyst)
View(tax_table(pscyst))
#removed ~4000 taxa

psroot <- prune_taxa(taxa_sums(psroot) > 0, psroot)
View(taxa_sums(psroot))
View(tax_table(psroot))
#removed ~4000 taxa



#Alpha diversity: Observed, Shannon & Simpson
#Need to split the shannon and simpson plots to view the box plots without the individual dots. 
#https://stackoverflow.com/questions/73544659/geom-points-are-not-placed-on-the-boxplot

#ps3 object- all samples
ps3alpha <- estimate_richness(ps3, measures = c("Observed", "Shannon", "Simpson"))
ps3alpha$Sample_Type <- sample_data(ps3)$Sample_Type
ps3alpha$Plot <- sample_data(ps3)$Plot
ps3alpha$Season <- sample_data(ps3)$Season
ps3alpha$SoilType <- sample_data(ps3)$SoilType

ggplot(ps3alpha, aes(x = Plot, y = Observed, color = Season)) + geom_boxplot()+ facet_grid(~Sample_Type, scales = "free", space = "free" )
#statistics
ps3alpha.sampletype <- aov(Observed ~Sample_Type, ps3alpha)
anova(ps3alpha.sampletype)
ps3alpha.season <- aov(Observed ~Season, ps3alpha)
anova(ps3alpha.season)
ps3alpha.plot <- aov(Observed ~Plot, ps3alpha)
anova(ps3alpha.plot)

ggplot(ps3alpha, aes(x = Plot, y = Shannon, color = Season)) + geom_boxplot()+ facet_grid(~Sample_Type, scales = "free", space = "free" )
#statistics
ps3alpha.sampletype <- aov(Shannon ~Sample_Type, ps3alpha)
anova(ps3alpha.sampletype)
ps3alpha.season <- aov(Shannon ~Season, ps3alpha)
anova(ps3alpha.season)
ps3alpha.plot <- aov(Shannon ~Plot, ps3alpha)
anova(ps3alpha.plot)

ggplot(ps3alpha, aes(x = Plot, y = Simpson, color = Season)) + geom_boxplot()+ facet_grid(~Sample_Type, scales = "free", space = "free" )
#statistics
ps3alpha.sampletype <- aov(Simpson ~Sample_Type, ps3alpha)
anova(ps3alpha.sampletype)
ps3alpha.season <- aov(Simpson ~Season, ps3alpha)
anova(ps3alpha.season)
ps3alpha.plot <- aov(Simpson ~Plot, ps3alpha)
anova(ps3alpha.plot)


#Do a dunn's test when the anova is significant (aka the means are not equal)
library(dunn.test)
dunn.test(ps3alpha$Observed, g = ps3alpha$Sample_Type, method = "bonferroni")

dunn.test(ps3alpha$Shannon, g = ps3alpha$Sample_Type, method = "bonferroni")

dunn.test(ps3alpha$Simpson, g = ps3alpha$Sample_Type, method = "bonferroni")


#plot_richness(pssoil, x="Plot", measures = c("Shannon", "Simpson"), color = "SoilType", shape = "Season") 
soilalpha <- estimate_richness(pssoil, measures = c("Observed","Shannon", "Simpson"))
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
rootalpha.plot <- aov(Simpson ~Plot, rootalpha)
anova(rootalpha.plot)


#Beta Diversity. Using Unifrac

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


#Beta diversity- Unifrac Soil 
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

#beta dispersion
pssoiltu.bdisper <- betadisper(pssoiltunifrac_dist, pssoil.t.tree.data$SoilType, "centroid", bias.adjust = TRUE)
anova(pssoiltu.bdisper)
permutest(pssoiltu.bdisper, pairwise = TRUE)



#Beta diversity- Unifrac mid season soil ---- 
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
pairwiseAdonis::pairwise.adonis(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
pssoiltmidu.bdisper <- betadisper(pssoil.midtunifrac_dist, pssoil.mid.t.tree.data$SoilType, "centroid", bias.adjust = TRUE)
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

#beta dispersion
psrhiztu.bdisper <- betadisper(psrhiztunifrac_dist, psrhiz.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psrhiztu.bdisper)
permutest(psrhiztu.bdisper, pairwise = TRUE)



#Beta diversity- Unifrac Bulk 
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
pairwiseAdonis::pairwise.adonis(psbulktunifrac_dist, psbulk.t.tree.data$Season, p.adjust.m = "holm", perm = 9999)

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
pairwiseAdonis::pairwise.adonis(psbulktunifrac_dist, psbulk.t.tree.data$Season, p.adjust.m = "holm", perm = 9999)

pairwiseAdonis::pairwise.adonis(psbulktunifrac_dist, psbulk.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#beta dispersion
psbulktu.bdisper <- betadisper(psbulktunifrac_dist, psbulk.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulktu.bdisper)
permutest(psbulktu.bdisper, pairwise = TRUE)



#Beta diversity- Unifrac Bulk Mid season samples 
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

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psbulkmidtunifrac_dist, psbulkmid.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

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

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psbulkmidtunifrac_dist, psbulkmid.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
psbulkmidtu.bdisper <- betadisper(psbulkmidtunifrac_dist, psbulkmid.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psbulkmidtu.bdisper)
permutest(psbulkmidtu.bdisper, pairwise = TRUE)


#Beta diversity- Unifrac Bulk Fall season samples 
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


#Beta diversity- Unifrac Cyst 
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

#betadispersion
pscysttu.bdisper <- betadisper(pscysttunifrac_dist, pscyst.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(pscysttu.bdisper)
permutest(pscysttu.bdisper, pairwise = TRUE)


#Beta diversity- Unifrac Cysts mid season samples 
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


#Beta diversity- Unifrac Cysts fall season samples 
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


#Beta diversity- Unifrac Roots 
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

#permanova posthoc
pairwiseAdonis::pairwise.adonis(psroottunifrac_dist, psroot.t.tree.data$Plot, p.adjust.m = "holm", perm = 9999)

#betadispersion
psroottu.bdisper <- betadisper(psroottunifrac_dist, psroot.t.tree.data$Plot, "centroid", bias.adjust = TRUE)
anova(psroottu.bdisper)
permutest(psroottu.bdisper, pairwise = TRUE)


## Make phyloseq objects with relative abundance values 
psrhiz.rab <- transform_sample_counts(psrhiz, function(x) x / sum(x))

psbulk.rab <- transform_sample_counts(psbulk, function(x) x / sum(x))

psbulkmid.rab <- transform_sample_counts(psbulkmid, function(x) x / sum(x))

psbulkfall.rab <- transform_sample_counts(psbulkfall, function(x) x / sum(x))

pscyst.rab <- transform_sample_counts(pscyst, function(x) x/sum(x))

pscystmid.rab <- transform_sample_counts(pscystmid, function(x) x/sum(x))

pscystfall.rab <- transform_sample_counts(pscystfall, function(x) x/sum(x))

psroot.rab <- transform_sample_counts(psroot, function(x) x/sum(x))


### Split ITS bulk soil into individual objects 
psbulkSA.ITS <- subset_samples(psbulk, Plot == "Sa")

psbulkS3.ITS <- subset_samples(psbulk, Plot == "S3")

psbulkSS.ITS <- subset_samples(psbulk, Plot == "Ss")

psbulkSR.ITS <- subset_samples(psbulk, Plot == "Sr")

##### split into mid season
psbulkSAmid.ITS <- subset_samples(psbulkSA.ITS, Season == "Mid")

psbulkS3mid.ITS <- subset_samples(psbulkS3.ITS, Season == "Mid")

psbulkSSmid.ITS <- subset_samples(psbulkSS.ITS, Season == "Mid")

psbulkSRmid.ITS <- subset_samples(psbulkSR.ITS, Season == "Mid")

#### split into fall season
psbulkSAfall.ITS <- subset_samples(psbulkSA.ITS, Season == "Fall")

psbulkS3fall.ITS <- subset_samples(psbulkS3.ITS, Season == "Fall")

psbulkSSfall.ITS <- subset_samples(psbulkSS.ITS, Season == "Fall")

psbulkSRfall.ITS <- subset_samples(psbulkSR.ITS, Season == "Fall")

##### remove taxa that have 0 sums.
psbulkSAmid.ITS <- prune_taxa(taxa_sums(psbulkSAmid.ITS) > 0, psbulkSAmid.ITS)

psbulkS3mid.ITS <- prune_taxa(taxa_sums(psbulkS3mid.ITS) > 0, psbulkS3mid.ITS)

psbulkSRmid.ITS <- prune_taxa(taxa_sums(psbulkSRmid.ITS) > 0, psbulkSRmid.ITS)

psbulkSSmid.ITS <- prune_taxa(taxa_sums(psbulkSSmid.ITS) > 0, psbulkSSmid.ITS)

psbulkSAfall.ITS <- prune_taxa(taxa_sums(psbulkSAfall.ITS) > 0, psbulkSAfall.ITS)

psbulkS3fall.ITS <- prune_taxa(taxa_sums(psbulkS3fall.ITS) > 0, psbulkS3fall.ITS)

psbulkSRfall.ITS <- prune_taxa(taxa_sums(psbulkSRfall.ITS) > 0, psbulkSRfall.ITS)

psbulkSSfall.ITS <- prune_taxa(taxa_sums(psbulkSSfall.ITS) > 0, psbulkSSfall.ITS)


### Split ITS cysts into individual objects 
pscystSA.ITS <- subset_samples(pscyst, Plot == "Sa")

pscystS3.ITS <- subset_samples(pscyst, Plot == "S3")

pscystSS.ITS <- subset_samples(pscyst, Plot == "Ss")

pscystSR.ITS <- subset_samples(pscyst, Plot == "Sr")

#### split into mid season
pscystSAmid.ITS <- subset_samples(pscystSA.ITS, Season == "Mid")

pscystS3mid.ITS <- subset_samples(pscystS3.ITS, Season == "Mid")

pscystSSmid.ITS <- subset_samples(pscystSS.ITS, Season == "Mid")

pscystSRmid.ITS <- subset_samples(pscystSR.ITS, Season == "Mid")

#### split into fall season
pscystSAfall.ITS <- subset_samples(pscystSA.ITS, Season == "Fall")

pscystS3fall.ITS <- subset_samples(pscystS3.ITS, Season == "Fall")

pscystSSfall.ITS <- subset_samples(pscystSS.ITS, Season == "Fall")

pscystSRfall.ITS <- subset_samples(pscystSR.ITS, Season == "Fall")

#### remove taxa that have 0 sums.
pscystSAmid.ITS <- prune_taxa(taxa_sums(pscystSAmid.ITS) > 0, pscystSAmid.ITS)

pscystS3mid.ITS <- prune_taxa(taxa_sums(pscystS3mid.ITS) > 0, pscystS3mid.ITS)

pscystSRmid.ITS <- prune_taxa(taxa_sums(pscystSRmid.ITS) > 0, pscystSRmid.ITS)

pscystSSmid.ITS <- prune_taxa(taxa_sums(pscystSSmid.ITS) > 0, pscystSSmid.ITS)


pscystSAfall.ITS <- prune_taxa(taxa_sums(pscystSAfall.ITS) > 0, pscystSAfall.ITS)

pscystS3fall.ITS <- prune_taxa(taxa_sums(pscystS3fall.ITS) > 0, pscystS3fall.ITS)

pscystSRfall.ITS <- prune_taxa(taxa_sums(pscystSRfall.ITS) > 0, pscystSRfall.ITS)

pscystSSfall.ITS <- prune_taxa(taxa_sums(pscystSSfall.ITS) > 0, pscystSSfall.ITS)

### Split ITS Rhizosphere into individual objects 
psrhizSA.ITS <- subset_samples(psrhiz, Plot == "Sa")

psrhizS3.ITS <- subset_samples(psrhiz, Plot == "S3")

psrhizSS.ITS <- subset_samples(psrhiz, Plot == "Ss")

psrhizSR.ITS <- subset_samples(psrhiz, Plot == "Sr")

#### remove taxa that have 0 sums.
psrhizSA.ITS <- prune_taxa(taxa_sums(psrhizSA.ITS) > 0, psrhizSA.ITS)

psrhizS3.ITS <- prune_taxa(taxa_sums(psrhizS3.ITS) > 0, psrhizS3.ITS)

psrhizSR.ITS <- prune_taxa(taxa_sums(psrhizSR.ITS) > 0, psrhizSR.ITS)

psrhizSS.ITS <- prune_taxa(taxa_sums(psrhizSS.ITS) > 0, psrhizSS.ITS)

### Split ITS Root into individual objects 
psrootSA.ITS <- subset_samples(psroot, Plot == "Sa")

psrootS3.ITS <- subset_samples(psroot, Plot == "S3")

psrootSS.ITS <- subset_samples(psroot, Plot == "Ss")

psrootSR.ITS <- subset_samples(psroot, Plot == "Sr")

#### remove taxa that have 0 sums.
psrootSA.ITS <- prune_taxa(taxa_sums(psrootSA.ITS) > 0, psrootSA.ITS)

psrootS3.ITS <- prune_taxa(taxa_sums(psrootS3.ITS) > 0, psrootS3.ITS)

psrootSR.ITS <- prune_taxa(taxa_sums(psrootSR.ITS) > 0, psrootSR.ITS)

psrootSS.ITS <- prune_taxa(taxa_sums(psrootSS.ITS) > 0, psrootSS.ITS)

## Cysts ITS plots- taxa with 1% abundance or higher- split by season 
#To make taxa be listed from high abundance to low abundance

#in ggplot command instead of fill = Genus

#replace with: factor(Genus, levels = c(setdiff(Genus, "Other"), "Other"))

### S3 mid season
pscystS3mid.ITS.rab <- transform_sample_counts(pscystS3mid.ITS, function(x) x/sum(x))

filterpscystS3mid.ITS.rab <- phyloseq::genefilter_sample(pscystS3mid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystS3mid.ITScore.rab <- prune_taxa(filterpscystS3mid.ITS.rab, pscystS3mid.ITS.rab)

cystS3mid.ITScore.genus.rab <- tax_glom(cystS3mid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#6 NA are being combined. Need to make sure they stay separate

cystS3mid.ITScore.genus.rab <-  renameTaxa(cystS3mid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystS3mid.ITScore.genus.rab.df <- psmelt(cystS3mid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystS3mid.ITScore.genus.rab.df$Genus <- as.character(cystS3mid.ITScore.genus.rab.df$Genus)

ggplot(cystS3mid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA mid season
pscystSAmid.ITS.rab <- transform_sample_counts(pscystSAmid.ITS, function(x) x/sum(x))

filterpscystSAmid.ITS.rab <- phyloseq::genefilter_sample(pscystSAmid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSAmid.ITScore.rab <- prune_taxa(filterpscystSAmid.ITS.rab, pscystSAmid.ITS.rab)

cystSAmid.ITScore.genus.rab <- tax_glom(cystSAmid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#5 NA are being combined. Need to make sure they stay separate

cystSAmid.ITScore.genus.rab <-  renameTaxa(cystSAmid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSAmid.ITScore.genus.rab.df <- psmelt(cystSAmid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSAmid.ITScore.genus.rab.df$Genus <- as.character(cystSAmid.ITScore.genus.rab.df$Genus)

ggplot(cystSAmid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR mid season
pscystSRmid.ITS.rab <- transform_sample_counts(pscystSRmid.ITS, function(x) x/sum(x))

filterpscystSRmid.ITS.rab <- phyloseq::genefilter_sample(pscystSRmid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSRmid.ITScore.rab <- prune_taxa(filterpscystSRmid.ITS.rab, pscystSRmid.ITS.rab)

cystSRmid.ITScore.genus.rab <- tax_glom(cystSRmid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#3 NA are being combined. Need to make sure they stay separate

cystSRmid.ITScore.genus.rab <-  renameTaxa(cystSRmid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe

cystSRmid.ITScore.genus.rab.df <- psmelt(cystSRmid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSRmid.ITScore.genus.rab.df$Genus <- as.character(cystSRmid.ITScore.genus.rab.df$Genus)

ggplot(cystSRmid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS mid season
pscystSSmid.ITS.rab <- transform_sample_counts(pscystSSmid.ITS, function(x) x/sum(x))

filterpscystSSmid.ITS.rab <- phyloseq::genefilter_sample(pscystSSmid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSSmid.ITScore.rab <- prune_taxa(filterpscystSSmid.ITS.rab, pscystSSmid.ITS.rab)

cystSSmid.ITScore.genus.rab <- tax_glom(cystSSmid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#8 NA are being combined. Need to make sure they stay separate

cystSSmid.ITScore.genus.rab <-  renameTaxa(cystSSmid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe

cystSSmid.ITScore.genus.rab.df <- psmelt(cystSSmid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSSmid.ITScore.genus.rab.df$Genus <- as.character(cystSSmid.ITScore.genus.rab.df$Genus)

ggplot(cystSSmid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### S3 fall season
pscystS3fall.ITS.rab <- transform_sample_counts(pscystS3fall.ITS, function(x) x/sum(x))

filterpscystS3fall.ITS.rab <- phyloseq::genefilter_sample(pscystS3fall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystS3fall.ITScore.rab <- prune_taxa(filterpscystS3fall.ITS.rab, pscystS3fall.ITS.rab)

cystS3fall.ITScore.genus.rab <- tax_glom(cystS3fall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#2 NA was unlabeled.

cystS3fall.ITScore.genus.rab <-  renameTaxa(cystS3fall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystS3fall.ITScore.genus.rab.df <- psmelt(cystS3fall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystS3fall.ITScore.genus.rab.df$Genus <- as.character(cystS3fall.ITScore.genus.rab.df$Genus)

ggplot(cystS3fall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA fall season
pscystSAfall.ITS.rab <- transform_sample_counts(pscystSAfall.ITS, function(x) x/sum(x))

filterpscystSAfall.ITS.rab <- phyloseq::genefilter_sample(pscystSAfall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSAfall.ITScore.rab <- prune_taxa(filterpscystSAfall.ITS.rab, pscystSAfall.ITS.rab)

cystSAfall.ITScore.genus.rab <- tax_glom(cystSAfall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#5 NA was unlabeled.

cystSAfall.ITScore.genus.rab <-  renameTaxa(cystSAfall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSAfall.ITScore.genus.rab.df <- psmelt(cystSAfall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSAfall.ITScore.genus.rab.df$Genus <- as.character(cystSAfall.ITScore.genus.rab.df$Genus)

ggplot(cystSAfall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR fall season
pscystSRfall.ITS.rab <- transform_sample_counts(pscystSRfall.ITS, function(x) x/sum(x))

filterpscystSRfall.ITS.rab <- phyloseq::genefilter_sample(pscystSRfall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSRfall.ITScore.rab <- prune_taxa(filterpscystSRfall.ITS.rab, pscystSRfall.ITS.rab)

cystSRfall.ITScore.genus.rab <- tax_glom(cystSRfall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#6 NA are being combined. Need to make sure they stay separate

cystSRfall.ITScore.genus.rab <-  renameTaxa(cystSRfall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSRfall.ITScore.genus.rab.df <- psmelt(cystSRfall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSRfall.ITScore.genus.rab.df$Genus <- as.character(cystSRfall.ITScore.genus.rab.df$Genus)

ggplot(cystSRfall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS fall season
pscystSSfall.ITS.rab <- transform_sample_counts(pscystSSfall.ITS, function(x) x/sum(x))

filterpscystSSfall.ITS.rab <- phyloseq::genefilter_sample(pscystSSfall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSSfall.ITScore.rab <- prune_taxa(filterpscystSSfall.ITS.rab, pscystSSfall.ITS.rab)

cystSSfall.ITScore.genus.rab <- tax_glom(cystSSfall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#3 NA are being combined. Need to make sure they stay separate

cystSSfall.ITScore.genus.rab <-  renameTaxa(cystSSfall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSSfall.ITScore.genus.rab.df <- psmelt(cystSSfall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSSfall.ITScore.genus.rab.df$Genus <- as.character(cystSSfall.ITScore.genus.rab.df$Genus)

ggplot(cystSSfall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


## Merge the Cyst 1% abundance or higher dataframes top taxa together then plot heatmap 
allcystmid.ITScore.genus.df <- rbind(cystS3mid.ITScore.genus.rab.df, cystSAmid.ITScore.genus.rab.df, cystSRmid.ITScore.genus.rab.df, cystSSmid.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

write.csv(allcystmid.ITScore.genus.df, "allcystmid.ITScore.genus.csv")

#import back in

allcystmid.ITScore.genus.df <- read.csv("fix.allcystmid.ITScore.genus.csv")

#ggplot(allcystmid.ITScore.genus.df, aes(x = Sample, y = Genus, fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank()) + facet_grid(~Plot, scales = "free", space = "free") + scale_fill_gradient(low = "yellow", high = "blue")
#ggplot(allcystmid.ITScore.genus.df, aes(x = Sample, y = Genus, fill = Abundance)) + geom_tile(color = "black"
) + theme(axis.text.x = element_blank()) + facet_grid(~Plot, scales = "free", space = "free") + scale_fill_gradientn(colors = hcl.colors(10, "Oslo"))

#make legend max, 84%
ggplot(allcystmid.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = 0.42
                                   , limits=c(0.0, 0.84), breaks=seq(0.0, 0.84, by=0.14))

allcystfall.ITScore.genus.df <- rbind(cystS3fall.ITScore.genus.rab.df, cystSAfall.ITScore.genus.rab.df, cystSRfall.ITScore.genus.rab.df, cystSSfall.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

write.csv(allcystfall.ITScore.genus.df, "allcystfall.ITScore.genus.csv")

#import back in

allcystfall.ITScore.genus.df <- read.csv("fix.allcystfall.ITScore.genus.csv")

#max legend 95%

ggplot(allcystfall.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black"
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"),axis.text.y.left = element_text(size =11),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = 0.475
                                   , limits=c(0.0, 0.95), breaks=seq(0.0, 0.95 ,by=0.09))


#ggplot(allcystfall.ITScore.genus.df, aes(x = Sample, y = Genus, fill = Abundance)) + geom_tile(color = "black", ) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y.left = element_text(size =11)) + facet_grid(~Plot, scales = "free", space = "free") + scale_fill_gradient(low = "blue",  high = "orange")



## Cyst ITS plots- taxa with 1% abundance or higher- NOT split by season 
#pscystSA.ITS <- subset_samples(pscyst, Plot == "Sa")

#pscystS3.ITS <- subset_samples(pscyst, Plot == "S3")

#pscystSS.ITS <- subset_samples(pscyst, Plot == "Ss")

#pscystSR.ITS <- subset_samples(pscyst, Plot == "Sr")

### S3  
pscystS3.ITS.rab <- transform_sample_counts(pscystS3.ITS, function(x) x/sum(x))

filterpscystS3.ITS.rab <- phyloseq::genefilter_sample(pscystS3.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystS3.ITScore.rab <- prune_taxa(filterpscystS3.ITS.rab, pscystS3.ITS.rab)

cystS3.ITScore.genus.rab <- tax_glom(cystS3.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#6 NA are being combined. Need to make sure they stay separate

cystS3.ITScore.genus.rab <-  renameTaxa(cystS3.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystS3.ITScore.genus.rab.df <- psmelt(cystS3.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystS3.ITScore.genus.rab.df$Genus <- as.character(cystS3.ITScore.genus.rab.df$Genus)

ggplot(cystS3.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA  
pscystSA.ITS.rab <- transform_sample_counts(pscystSA.ITS, function(x) x/sum(x))

filterpscystSA.ITS.rab <- phyloseq::genefilter_sample(pscystSA.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSA.ITScore.rab <- prune_taxa(filterpscystSA.ITS.rab, pscystSA.ITS.rab)

cystSA.ITScore.genus.rab <- tax_glom(cystSA.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#9 NA was being combined. Need to make sure they stay separate

cystSA.ITScore.genus.rab <-  renameTaxa(cystSA.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSA.ITScore.genus.rab.df <- psmelt(cystSA.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSA.ITScore.genus.rab.df$Genus <- as.character(cystSA.ITScore.genus.rab.df$Genus)

ggplot(cystSA.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR  
pscystSR.ITS.rab <- transform_sample_counts(pscystSR.ITS, function(x) x/sum(x))

filterpscystSR.ITS.rab <- phyloseq::genefilter_sample(pscystSR.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSR.ITScore.rab <- prune_taxa(filterpscystSR.ITS.rab, pscystSR.ITS.rab)

cystSR.ITScore.genus.rab <- tax_glom(cystSR.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#7 NA are being combined. Need to make sure they stay separate

cystSR.ITScore.genus.rab <-  renameTaxa(cystSR.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSR.ITScore.genus.rab.df <- psmelt(cystSR.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSR.ITScore.genus.rab.df$Genus <- as.character(cystSR.ITScore.genus.rab.df$Genus)

ggplot(cystSR.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS  
pscystSS.ITS.rab <- transform_sample_counts(pscystSS.ITS, function(x) x/sum(x))

filterpscystSS.ITS.rab <- phyloseq::genefilter_sample(pscystSS.ITS.rab, filterfun_sample(function(x) x / sum(x) > .01))

cystSS.ITScore.rab <- prune_taxa(filterpscystSS.ITS.rab, pscystSS.ITS.rab)

cystSS.ITScore.genus.rab <- tax_glom(cystSS.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#9 NA are being combined. Need to make sure they stay separate

cystSS.ITScore.genus.rab <-  renameTaxa(cystSS.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

cystSS.ITScore.genus.rab.df <- psmelt(cystSS.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

cystSS.ITScore.genus.rab.df$Genus <- as.character(cystSS.ITScore.genus.rab.df$Genus)

ggplot(cystSS.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


## Merge Cyst 1% abundance or higher for Heatmap - Not split by season 
allcysts.ITScore.genus.df <- rbind(cystS3.ITScore.genus.rab.df, cystSA.ITScore.genus.rab.df, cystSR.ITScore.genus.rab.df, cystSS.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

write.csv(allcysts.ITScore.genus.df, "allcysts.ITScore.genus.csv")

#import back in, add prevalence & total abundance 

allcysts.ITScore.genus.df <- read.csv("fix.allcysts.ITScore.genus.csv")

allcysts.ITScore.genus.df$Season <- factor(allcysts.ITScore.genus.prev.df$Season, levels = c("Mid", "Fall"))

#WINNER max 94.6% 

#need to change the values depending on how you want to view it. 

ggplot(allcysts.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.20, 0.40, 0.60, 0.80, 0.946), values = c("0", "0.2", "1"))


## Bulk ITS- treatments- not split by season - 2% or higher
#psbulkSA.ITS <- subset_samples(psbulk, Plot == "Sa")

#psbulkS3.ITS <- subset_samples(psbulk, Plot == "S3")

#psbulkSS.ITS <- subset_samples(psbulk, Plot == "Ss")

#psbulkSR.ITS <- subset_samples(psbulk, Plot == "Sr")


### S3  
psbulkS3.ITS.rab <- transform_sample_counts(psbulkS3.ITS, function(x) x/sum(x))

filterpsbulkS3.ITS.rab <- phyloseq::genefilter_sample(psbulkS3.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkS3.ITScore.rab <- prune_taxa(filterpsbulkS3.ITS.rab, psbulkS3.ITS.rab)

bulkS3.ITScore.genus.rab <- tax_glom(bulkS3.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#5 NA are being combined. Need to make sure they stay separate

bulkS3.ITScore.genus.rab <-  renameTaxa(bulkS3.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkS3.ITScore.genus.rab.df <- psmelt(bulkS3.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkS3.ITScore.genus.rab.df$Genus <- as.character(bulkS3.ITScore.genus.rab.df$Genus)

ggplot(bulkS3.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA  
psbulkSA.ITS.rab <- transform_sample_counts(psbulkSA.ITS, function(x) x/sum(x))

filterpsbulkSA.ITS.rab <- phyloseq::genefilter_sample(psbulkSA.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSA.ITScore.rab <- prune_taxa(filterpsbulkSA.ITS.rab, psbulkSA.ITS.rab)

bulkSA.ITScore.genus.rab <- tax_glom(bulkSA.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#4 NA was being combined. Need to make sure they stay separate

bulkSA.ITScore.genus.rab <-  renameTaxa(bulkSA.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkSA.ITScore.genus.rab.df <- psmelt(bulkSA.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSA.ITScore.genus.rab.df$Genus <- as.character(bulkSA.ITScore.genus.rab.df$Genus)

ggplot(bulkSA.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR  
psbulkSR.ITS.rab <- transform_sample_counts(psbulkSR.ITS, function(x) x/sum(x))

filterpsbulkSR.ITS.rab <- phyloseq::genefilter_sample(psbulkSR.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSR.ITScore.rab <- prune_taxa(filterpsbulkSR.ITS.rab, psbulkSR.ITS.rab)

bulkSR.ITScore.genus.rab <- tax_glom(bulkSR.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#6 NA are being combined. Need to make sure they stay separate

bulkSR.ITScore.genus.rab <-  renameTaxa(bulkSR.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkSR.ITScore.genus.rab.df <- psmelt(bulkSR.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSR.ITScore.genus.rab.df$Genus <- as.character(bulkSR.ITScore.genus.rab.df$Genus)

ggplot(bulkSR.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS  
psbulkSS.ITS.rab <- transform_sample_counts(psbulkSS.ITS, function(x) x/sum(x))

filterpsbulkSS.ITS.rab <- phyloseq::genefilter_sample(psbulkSS.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSS.ITScore.rab <- prune_taxa(filterpsbulkSS.ITS.rab, psbulkSS.ITS.rab)

bulkSS.ITScore.genus.rab <- tax_glom(bulkSS.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#5 NA are being combined. Need to make sure they stay separate

bulkSS.ITScore.genus.rab <-  renameTaxa(bulkSS.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkSS.ITScore.genus.rab.df <- psmelt(bulkSS.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSS.ITScore.genus.rab.df$Genus <- as.character(bulkSS.ITScore.genus.rab.df$Genus)

ggplot(bulkSS.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


## Merge Bulk ITS dataframes together for Heatmap of 2% or higher- not split by season 
allbulk.ITScore.genus.df <- rbind(bulkS3.ITScore.genus.rab.df, bulkSA.ITScore.genus.rab.df, bulkSR.ITScore.genus.rab.df, bulkSS.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

#added prevlence and total abundance too

write.csv(allbulk.ITScore.genus.df, "allbulk.ITScore.genus.csv")

#import back in

allbulk.ITScore.genus.df <- read.csv("fix.allbulk.ITScore.genus.csv")

#make mid season first

allbulk.ITScore.genus.df$Season <- factor(allbulk.ITScore.genus.df$Season, levels = c("Mid", "Fall"))

#WINNER 49.2%

#need to change the values depending on how you want to view it. 

ggplot(allbulk.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.492), values = c("0", "0.2", "1"))


## Bulk Soil ITS plots- taxa with 2% abundance or higher- split by season 
#To make taxa be listed from high abundance to low abundance

#in ggplot command instead of fill = Genus

#replace with: factor(Genus, levels = c(setdiff(Genus, "Other"), "Other"))

### S3 mid season
psbulkS3mid.ITS.rab <- transform_sample_counts(psbulkS3mid.ITS, function(x) x/sum(x))

filterpsbulkS3mid.ITS.rab <- phyloseq::genefilter_sample(psbulkS3mid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkS3mid.ITScore.rab <- prune_taxa(filterpsbulkS3mid.ITS.rab, psbulkS3mid.ITS.rab)

bulkS3mid.ITScore.genus.rab <- tax_glom(bulkS3mid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#4 NA are being combined. Need to make sure they stay separate

bulkS3mid.ITScore.genus.rab <-  renameTaxa(bulkS3mid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkS3mid.ITScore.genus.rab.df <- psmelt(bulkS3mid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkS3mid.ITScore.genus.rab.df$Genus <- as.character(bulkS3mid.ITScore.genus.rab.df$Genus)

ggplot(bulkS3mid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA mid season
psbulkSAmid.ITS.rab <- transform_sample_counts(psbulkSAmid.ITS, function(x) x/sum(x))

filterpsbulkSAmid.ITS.rab <- phyloseq::genefilter_sample(psbulkSAmid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSAmid.ITScore.rab <- prune_taxa(filterpsbulkSAmid.ITS.rab, psbulkSAmid.ITS.rab)

bulkSAmid.ITScore.genus.rab <- tax_glom(bulkSAmid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#4 NAs are being combined. Need to make sure they stay separate

bulkSAmid.ITScore.genus.rab <-  renameTaxa(bulkSAmid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkSAmid.ITScore.genus.rab.df <- psmelt(bulkSAmid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSAmid.ITScore.genus.rab.df$Genus <- as.character(bulkSAmid.ITScore.genus.rab.df$Genus)

ggplot(bulkSAmid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR mid season
psbulkSRmid.ITS.rab <- transform_sample_counts(psbulkSRmid.ITS, function(x) x/sum(x))

filterpsbulkSRmid.ITS.rab <- phyloseq::genefilter_sample(psbulkSRmid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSRmid.ITScore.rab <- prune_taxa(filterpsbulkSRmid.ITS.rab, psbulkSRmid.ITS.rab)

bulkSRmid.ITScore.genus.rab <- tax_glom(bulkSRmid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#4 NA are being combined. Need to make sure they stay separate

bulkSRmid.ITScore.genus.rab <-  renameTaxa(bulkSRmid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkSRmid.ITScore.genus.rab.df <- psmelt(bulkSRmid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSRmid.ITScore.genus.rab.df$Genus <- as.character(bulkSRmid.ITScore.genus.rab.df$Genus)

ggplot(bulkSRmid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS mid season
psbulkSSmid.ITS.rab <- transform_sample_counts(psbulkSSmid.ITS, function(x) x/sum(x))

filterpsbulkSSmid.ITS.rab <- phyloseq::genefilter_sample(psbulkSSmid.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSSmid.ITScore.rab <- prune_taxa(filterpsbulkSSmid.ITS.rab, psbulkSSmid.ITS.rab)

bulkSSmid.ITScore.genus.rab <- tax_glom(bulkSSmid.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#3 NA are being combined. Need to make sure they stay separate

bulkSSmid.ITScore.genus.rab <-  renameTaxa(bulkSSmid.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                           numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe

bulkSSmid.ITScore.genus.rab.df <- psmelt(bulkSSmid.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSSmid.ITScore.genus.rab.df$Genus <- as.character(bulkSSmid.ITScore.genus.rab.df$Genus)

ggplot(bulkSSmid.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### S3 fall season
psbulkS3fall.ITS.rab <- transform_sample_counts(psbulkS3fall.ITS, function(x) x/sum(x))

filterpsbulkS3fall.ITS.rab <- phyloseq::genefilter_sample(psbulkS3fall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkS3fall.ITScore.rab <- prune_taxa(filterpsbulkS3fall.ITS.rab, psbulkS3fall.ITS.rab)

bulkS3fall.ITScore.genus.rab <- tax_glom(bulkS3fall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#2 NA was unlabeled.

bulkS3fall.ITScore.genus.rab <-  renameTaxa(bulkS3fall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe

bulkS3fall.ITScore.genus.rab.df <- psmelt(bulkS3fall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkS3fall.ITScore.genus.rab.df$Genus <- as.character(bulkS3fall.ITScore.genus.rab.df$Genus)

ggplot(bulkS3fall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA fall season
psbulkSAfall.ITS.rab <- transform_sample_counts(psbulkSAfall.ITS, function(x) x/sum(x))

filterpsbulkSAfall.ITS.rab <- phyloseq::genefilter_sample(psbulkSAfall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSAfall.ITScore.rab <- prune_taxa(filterpsbulkSAfall.ITS.rab, psbulkSAfall.ITS.rab)

bulkSAfall.ITScore.genus.rab <- tax_glom(bulkSAfall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#1 NA was unlabeled.

bulkSAfall.ITScore.genus.rab <-  renameTaxa(bulkSAfall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe

bulkSAfall.ITScore.genus.rab.df <- psmelt(bulkSAfall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSAfall.ITScore.genus.rab.df$Genus <- as.character(bulkSAfall.ITScore.genus.rab.df$Genus)

ggplot(bulkSAfall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR fall season
psbulkSRfall.ITS.rab <- transform_sample_counts(psbulkSRfall.ITS, function(x) x/sum(x))

filterpsbulkSRfall.ITS.rab <- phyloseq::genefilter_sample(psbulkSRfall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSRfall.ITScore.rab <- prune_taxa(filterpsbulkSRfall.ITS.rab, psbulkSRfall.ITS.rab)

bulkSRfall.ITScore.genus.rab <- tax_glom(bulkSRfall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#2 NA are being combined. Need to make sure they stay separate

bulkSRfall.ITScore.genus.rab <-  renameTaxa(bulkSRfall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

bulkSRfall.ITScore.genus.rab.df <- psmelt(bulkSRfall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSRfall.ITScore.genus.rab.df$Genus <- as.character(bulkSRfall.ITScore.genus.rab.df$Genus)

ggplot(bulkSRfall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS fall season
psbulkSSfall.ITS.rab <- transform_sample_counts(psbulkSSfall.ITS, function(x) x/sum(x))

filterpsbulkSSfall.ITS.rab <- phyloseq::genefilter_sample(psbulkSSfall.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

bulkSSfall.ITScore.rab <- prune_taxa(filterpsbulkSSfall.ITS.rab, psbulkSSfall.ITS.rab)

bulkSSfall.ITScore.genus.rab <- tax_glom(bulkSSfall.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#four NA are being combined. Need to make sure they stay separate

bulkSSfall.ITScore.genus.rab <-  renameTaxa(bulkSSfall.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                            numDupli = "Genus", numDupliPat ="<name><.><num>")
#make dataframe

bulkSSfall.ITScore.genus.rab.df <- psmelt(bulkSSfall.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

bulkSSfall.ITScore.genus.rab.df$Genus <- as.character(bulkSSfall.ITScore.genus.rab.df$Genus)

ggplot(bulkSSfall.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")



## Merge Bulk soil ITS dataframes together and heatmap of the 2% abundant- Split by season 
### Bulk Mid
allbulkmid.ITScore.genus.df <- rbind(bulkS3mid.ITScore.genus.rab.df, bulkSAmid.ITScore.genus.rab.df, bulkSRmid.ITScore.genus.rab.df, bulkSSmid.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

write.csv(allbulkmid.ITScore.genus.df, "allbulkmid.ITScore.genus.csv")

#import back in

allbulkmid.ITScore.genus.df <- read.csv("fix.allbulkmid.ITScore.genus.csv")

#label max legend as 50%

ggplot(allbulkmid.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(Abundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"),axis.text.y.left = element_text(size =10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "white", mid = "orange", high = "blue", midpoint = .25
                                   , limits=c(0.0, 0.50), breaks=seq(0.0, 0.50 ,by=0.1))


### Bulk fall 
allbulkfall.ITScore.genus.df <- rbind(bulkS3fall.ITScore.genus.rab.df, bulkSAfall.ITScore.genus.rab.df, bulkSRfall.ITScore.genus.rab.df, bulkSSfall.ITScore.genus.rab.df)
#fix names of family taxa to have the same name. 
write.csv(allbulkfall.ITScore.genus.df, "allbulkfall.ITScore.genus.csv")
#import back in
allbulkfall.ITScore.genus.df <- read.csv("fix.allbulkfall.ITScore.genus.csv")

ggplot(allbulkfall.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.20, 0.40, 0.60, 0.80, 0.946), values = c("0", "0.2", "1"))



## Rhiz Soil ITS plots- taxa with 2% abundance or higher 
#To make taxa be listed from high abundance to low abundance

#in ggplot command instead of fill = Genus

#replace with: factor(Genus, levels = c(setdiff(Genus, "Other"), "Other"))

### S3  season
psrhizS3.ITS.rab <- transform_sample_counts(psrhizS3.ITS, function(x) x/sum(x))

filterpsrhizS3.ITS.rab <- phyloseq::genefilter_sample(psrhizS3.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rhizS3.ITScore.rab <- prune_taxa(filterpsrhizS3.ITS.rab, psrhizS3.ITS.rab)

rhizS3.ITScore.genus.rab <- tax_glom(rhizS3.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#3 NA are being combined. Need to make sure they stay separate

rhizS3.ITScore.genus.rab <-  renameTaxa(rhizS3.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rhizS3.ITScore.genus.rab.df <- psmelt(rhizS3.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rhizS3.ITScore.genus.rab.df$Genus <- as.character(rhizS3.ITScore.genus.rab.df$Genus)

ggplot(rhizS3.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SA  season
psrhizSA.ITS.rab <- transform_sample_counts(psrhizSA.ITS, function(x) x/sum(x))

filterpsrhizSA.ITS.rab <- phyloseq::genefilter_sample(psrhizSA.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rhizSA.ITScore.rab <- prune_taxa(filterpsrhizSA.ITS.rab, psrhizSA.ITS.rab)

rhizSA.ITScore.genus.rab <- tax_glom(rhizSA.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#5 NA are being combined. Need to make sure they stay separate

rhizSA.ITScore.genus.rab <-  renameTaxa(rhizSA.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rhizSA.ITScore.genus.rab.df <- psmelt(rhizSA.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rhizSA.ITScore.genus.rab.df$Genus <- as.character(rhizSA.ITScore.genus.rab.df$Genus)

ggplot(rhizSA.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR  season
psrhizSR.ITS.rab <- transform_sample_counts(psrhizSR.ITS, function(x) x/sum(x))

filterpsrhizSR.ITS.rab <- phyloseq::genefilter_sample(psrhizSR.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rhizSR.ITScore.rab <- prune_taxa(filterpsrhizSR.ITS.rab, psrhizSR.ITS.rab)

rhizSR.ITScore.genus.rab <- tax_glom(rhizSR.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#3 NA are being combined. Need to make sure they stay separate

rhizSR.ITScore.genus.rab <-  renameTaxa(rhizSR.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rhizSR.ITScore.genus.rab.df <- psmelt(rhizSR.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rhizSR.ITScore.genus.rab.df$Genus <- as.character(rhizSR.ITScore.genus.rab.df$Genus)

ggplot(rhizSR.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS  season
psrhizSS.ITS.rab <- transform_sample_counts(psrhizSS.ITS, function(x) x/sum(x))

filterpsrhizSS.ITS.rab <- phyloseq::genefilter_sample(psrhizSS.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rhizSS.ITScore.rab <- prune_taxa(filterpsrhizSS.ITS.rab, psrhizSS.ITS.rab)

rhizSS.ITScore.genus.rab <- tax_glom(rhizSS.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#2 NA are being combined. Need to make sure they stay separate

rhizSS.ITScore.genus.rab <-  renameTaxa(rhizSS.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rhizSS.ITScore.genus.rab.df <- psmelt(rhizSS.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rhizSS.ITScore.genus.rab.df$Genus <- as.character(rhizSS.ITScore.genus.rab.df$Genus)

ggplot(rhizSS.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


## Merge Rhizosphere ITS dataframes together and make a heatmap of the 2% abundant or higher taxa ----
allrhiz.ITScore.genus.df <- rbind(rhizS3.ITScore.genus.rab.df, rhizSA.ITScore.genus.rab.df, rhizSR.ITScore.genus.rab.df, rhizSS.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

write.csv(allrhiz.ITScore.genus.df, "allrhiz.ITScore.genus.csv")

#import back in

allrhiz.ITScore.genus.df <- read.csv("fix.allrhiz.ITScore.genus.csv")

#winner 47.1%

ggplot(allrhiz.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.471), values = c("0", "0.2", "1"))

## Roots ITS taxa with 2% abundance or higher. 
### S3  season
psrootS3.ITS.rab <- transform_sample_counts(psrootS3.ITS, function(x) x/sum(x))

filterpsrootS3.ITS.rab <- phyloseq::genefilter_sample(psrootS3.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rootS3.ITScore.rab <- prune_taxa(filterpsrootS3.ITS.rab, psrootS3.ITS.rab)

rootS3.ITScore.genus.rab <- tax_glom(rootS3.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#6 NA are being combined. Need to make sure they stay separate

rootS3.ITScore.genus.rab <-  renameTaxa(rootS3.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rootS3.ITScore.genus.rab.df <- psmelt(rootS3.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rootS3.ITScore.genus.rab.df$Genus <- as.character(rootS3.ITScore.genus.rab.df$Genus)

ggplot(rootS3.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

### SA  season
psrootSA.ITS.rab <- transform_sample_counts(psrootSA.ITS, function(x) x/sum(x))

filterpsrootSA.ITS.rab <- phyloseq::genefilter_sample(psrootSA.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rootSA.ITScore.rab <- prune_taxa(filterpsrootSA.ITS.rab, psrootSA.ITS.rab)

rootSA.ITScore.genus.rab <- tax_glom(rootSA.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#8 NA are being combined. Need to make sure they stay separate

rootSA.ITScore.genus.rab <-  renameTaxa(rootSA.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rootSA.ITScore.genus.rab.df <- psmelt(rootSA.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rootSA.ITScore.genus.rab.df$Genus <- as.character(rootSA.ITScore.genus.rab.df$Genus)

ggplot(rootSA.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SR  season
psrootSR.ITS.rab <- transform_sample_counts(psrootSR.ITS, function(x) x/sum(x))

filterpsrootSR.ITS.rab <- phyloseq::genefilter_sample(psrootSR.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rootSR.ITScore.rab <- prune_taxa(filterpsrootSR.ITS.rab, psrootSR.ITS.rab)

rootSR.ITScore.genus.rab <- tax_glom(rootSR.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#4 NA are being combined. Need to make sure they stay separate

rootSR.ITScore.genus.rab <-  renameTaxa(rootSR.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rootSR.ITScore.genus.rab.df <- psmelt(rootSR.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rootSR.ITScore.genus.rab.df$Genus <- as.character(rootSR.ITScore.genus.rab.df$Genus)

ggplot(rootSR.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


### SS  season
psrootSS.ITS.rab <- transform_sample_counts(psrootSS.ITS, function(x) x/sum(x))

filterpsrootSS.ITS.rab <- phyloseq::genefilter_sample(psrootSS.ITS.rab, filterfun_sample(function(x) x / sum(x) > .02))

rootSS.ITScore.rab <- prune_taxa(filterpsrootSS.ITS.rab, psrootSS.ITS.rab)

rootSS.ITScore.genus.rab <- tax_glom(rootSS.ITScore.rab, taxrank = "Genus", NArm = FALSE)

#7 NA are being combined. Need to make sure they stay separate

rootSS.ITScore.genus.rab <-  renameTaxa(rootSS.ITScore.genus.rab, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                        numDupli = "Genus", numDupliPat ="<name><.><num>")

#make dataframe

rootSS.ITScore.genus.rab.df <- psmelt(rootSS.ITScore.genus.rab)

#convert Genus to a character vector from a factor because R

rootSS.ITScore.genus.rab.df$Genus <- as.character(rootSS.ITScore.genus.rab.df$Genus)

ggplot(rootSS.ITScore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")


## Merge root dataframes for heatmap of 2% or higher.  
allroot.ITScore.genus.df <- rbind(rootS3.ITScore.genus.rab.df, rootSA.ITScore.genus.rab.df, rootSR.ITScore.genus.rab.df, rootSS.ITScore.genus.rab.df)

#fix names of family taxa to have the same name. 

write.csv(allroot.ITScore.genus.df, "allroot.ITScore.genus.csv")

#import back in

allroot.ITScore.genus.df <- read.csv("fix.allroot.ITScore.genus.csv")

#winner 82.4%

ggplot(allroot.ITScore.genus.df, aes(x = Sample, y= reorder(Genus, as.integer(factor(TotalAbundance))), fill = Abundance)) + geom_tile(color = "black", 
) + theme(axis.text.x = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 10),
          strip.text.x.top = element_text(size = 12)) + facet_grid(~Plot*Season, scales = "free", space = "free"
          ) + scale_fill_gradientn(colors = c("white", "orange" ,"blue")
                                   , breaks = c(0.0, 0.2, 0.4, 0.6, 0.824), values = c("0", "0.2", "1"))


## Prevalence plots 
#export allcystmid.ITScore.genus.df as a csv.

#in excel, make columns for # of samples present in, Prevalence (# of samples present/total samples), Median abundance, Total abundance

#Total abundance wasn't working like it did for reordering in the 16S data

#I made the total abundance a number value so that it would be in the same order as the heatmap. 

#import back in. 

### cysts 
cyst.all.prevalence.df <- read.csv("CystAllITSPrevalence.csv")

cyst.all.prevalence.df$Genus <- as.character(cyst.all.prevalence.df$Genus)

cyst.all.prevalence.df$Genus <- factor(cyst.all.prevalence.df$Genus, levels=unique(cyst.all.prevalence.df$Genus))

#barplot

ggplot(cyst.all.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()


### bulk 
bulk.all.prevalence.df <- read.csv("BulkAllITSPrevalence.csv")

bulk.all.prevalence.df$Genus <- as.character(bulk.all.prevalence.df$Genus)

bulk.all.prevalence.df$Genus <- factor(bulk.all.prevalence.df$Genus, levels=unique(bulk.all.prevalence.df$Genus))

ggplot(bulk.all.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()

### rhizosphere
rhiz.prevalence.df <- read.csv("RhizITSPrevalence.csv")

rhiz.prevalence.df$Genus <- as.character(rhiz.prevalence.df$Genus)

rhiz.prevalence.df$Genus <- factor(rhiz.prevalence.df$Genus, levels=unique(rhiz.prevalence.df$Genus))

#ggplot(rhiz.prevalence.df, aes(x= Prevalence, y= reorder(Genus, as.integer(factor(HeatmapOrder)))
)) + geom_point(aes(size = Median)) + theme(axis.text.y.left = element_text(size = 12), axis.text.x.bottom = element_text(size = 10)) + ylab("Genus")

ggplot(rhiz.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()


### roots 
root.prevalence.df <- read.csv("RootITSPrevalence.csv")

root.prevalence.df$Genus <- as.character(root.prevalence.df$Genus)

root.prevalence.df$Genus <- factor(root.prevalence.df$Genus, levels=unique(root.prevalence.df$Genus))

#ggplot(root.prevalence.df, aes(x= Prevalence, y= reorder(Genus, as.integer(factor(HeatmapOrder)))
#)) + geom_point(aes(size = Median)) + theme(axis.text.y.left = element_text(size = 12), axis.text.x.bottom = element_text(size = 10)) + ylab("Genus")

ggplot(root.prevalence.df, aes(y= reorder(Genus, as.integer(factor(TotalAbundance))), x= Prevalence)) + geom_col()



## NetCoMi - cysts 
filterpscyst.ITS.rab <- phyloseq::genefilter_sample(pscyst.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

pscyst.ITS.rab.filt <- prune_taxa(filterpscyst.ITS.rab, pscyst.rab)

pscyst.ITS.rab.glom <- tax_glom(pscyst.ITS.rab.filt, taxrank = "Genus", NArm = FALSE)

pscyst.ITS.rab.rename <-  renameTaxa(pscyst.ITS.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_pscystglom_spiec <- netConstruct(pscyst.ITS.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_pscystglom_spiec <- netAnalyze(net_pscystglom_spiec)

summary(analyze_pscystglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_pscystglom_spiec, labelScale = TRUE, cexNodes = 1, labelFont = 3, rmSingles = TRUE, title1 = " Cysts glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi Cysts ITS Both Seasons 1percent.pdf" 


## NetCoMi- bulk soil 
filterpsbulk.ITS.rab <- phyloseq::genefilter_sample(psbulk.rab, filterfun_sample(function(x) x / sum(x) > 0.02))

psbulk.ITS.rab.filt <- prune_taxa(filterpsbulk.ITS.rab, psbulk.rab)

psbulk.ITS.rab.glom <- tax_glom(psbulk.ITS.rab.filt, taxrank = "Genus", NArm = FALSE)

psbulk.ITS.rab.rename <-  renameTaxa(psbulk.ITS.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psbulkglom_spiec <- netConstruct(psbulk.ITS.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psbulkglom_spiec <- netAnalyze(net_psbulkglom_spiec)

summary(analyze_psbulkglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psbulkglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Bulk glom >2% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi Bulk ITS Both Seasons 2percent.pdf"


## NetCoMi- Rhizosphere- Network is empty
filterpsrhiz.ITS.rab <- phyloseq::genefilter_sample(psrhiz.rab, filterfun_sample(function(x) x / sum(x) > 0.02))

psrhiz.ITS.rab.filt <- prune_taxa(filterpsrhiz.ITS.rab, psrhiz.rab)

psrhiz.ITS.rab.glom <- tax_glom(psrhiz.ITS.rab.filt, taxrank = "Genus", NArm = FALSE)

psrhiz.ITS.rab.rename <-  renameTaxa(psrhiz.ITS.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psrhizglom_spiec <- netConstruct(psrhiz.ITS.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psrhizglom_spiec <- netAnalyze(net_psrhizglom_spiec)

summary(analyze_psrhizglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psrhizglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Rhiz glom > 2% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)


## NetCoMi- roots 
filterpsroot.ITS.rab <- phyloseq::genefilter_sample(psroot.rab, filterfun_sample(function(x) x / sum(x) > 0.02))

psroot.ITS.rab.filt <- prune_taxa(filterpsroot.ITS.rab, psroot.rab)

psroot.ITS.rab.glom <- tax_glom(psroot.ITS.rab.filt, taxrank = "Genus", NArm = FALSE)

psroot.ITS.rab.rename <-  renameTaxa(psroot.ITS.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psrootglom_spiec <- netConstruct(psroot.ITS.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psrootglom_spiec <- netAnalyze(net_psrootglom_spiec)

summary(analyze_psrootglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psrootglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " Roots glom > 2% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi Roots ITS Both Seasons 2percent.pdf"



## NetCoMi- by treatment- all samples have the same cut off of 1%
ps3.rab <- transform_sample_counts(ps3, function(x) x/sum(x))

filterps3.rab <- phyloseq::genefilter_sample(ps3.rab, filterfun_sample(function(x) x / sum(x) > .01))

ps3.filt.rab <- prune_taxa(filterps3.rab, ps3.rab)

#split into treatments 

ps3.filt.S3 <- subset_samples(ps3.filt.rab, Plot == "S3")

ps3.filt.SA <- subset_samples(ps3.filt.rab, Plot == "Sa")

ps3.filt.SR <- subset_samples(ps3.filt.rab, Plot == "Sr")

ps3.filt.SS <- subset_samples(ps3.filt.rab, Plot == "Ss")


## NetCoMi by treatment-all samples 1% abund S3
ps3.filt.S3.glom <- tax_glom(ps3.filt.S3, taxrank = "Genus", NArm = FALSE)

ps3.filt.S3.rename <-  renameTaxa(ps3.filt.S3.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps3.S3_spiec <- netConstruct(ps3.filt.S3.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)
analyze_ps3.S3_spiec <- netAnalyze(net_ps3.S3_spiec)

summary(analyze_ps3.S3_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_ps3.S3_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "S3 ITS all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi S3 ITS all samples 1percent abund.pdf"


## NetCoMi by treatment-all samples 1% abund SA
ps3.filt.SA.glom <- tax_glom(ps3.filt.SA, taxrank = "Genus", NArm = FALSE)

ps3.filt.SA.rename <-  renameTaxa(ps3.filt.SA.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps3.SA_spiec <- netConstruct(ps3.filt.SA.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_ps3.SA_spiec <- netAnalyze(net_ps3.SA_spiec)

summary(analyze_ps3.SA_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_ps3.SA_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "SA ITS all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi SA ITS all samples 1percent abund.pdf"


## NetCoMi by treatment-all samples 1% abund SR
ps3.filt.SR.glom <- tax_glom(ps3.filt.SR, taxrank = "Genus", NArm = FALSE)

ps3.filt.SR.rename <-  renameTaxa(ps3.filt.SR.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps3.SR_spiec <- netConstruct(ps3.filt.SR.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_ps3.SR_spiec <- netAnalyze(net_ps3.SR_spiec)

summary(analyze_ps3.SR_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_ps3.SR_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "SR ITS all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi SR ITS all samples 1percent abund.pdf"


## NetCoMi by treatment-all samples 1% abund SS
ps3.filt.SS.glom <- tax_glom(ps3.filt.SS, taxrank = "Genus", NArm = FALSE)

ps3.filt.SS.rename <-  renameTaxa(ps3.filt.SS.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus", numDupliPat ="<name><.><num>")

net_ps3.SS_spiec <- netConstruct(ps3.filt.SS.rename, measure = "spieceasi", taxRank = "Genus", 
                                 measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                 sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_ps3.SS_spiec <- netAnalyze(net_ps3.SS_spiec)

summary(analyze_ps3.SS_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_ps3.SS_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = "SS ITS all 1%"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "NetCoMi SS ITS all samples 1percent abund.pdf"


## ANCOM-BC- entire dataset with thresholds from abundance plots 
library(ANCOMBC)

#Make a phyloseq object with the same taxa that are included in abundance thresholds. 

#Cysts = 1% or higher, Soils = 2% or higher, Roots = 2% or higher. 

pscyst.rab <- transform_sample_counts(pscyst, function(x) x/sum(x))
filter.pscyst.rab <- phyloseq::genefilter_sample(pscyst.rab, filterfun_sample(function(x) x / sum(x) > .01))
#pull from the non relative abundance bulk file
forancom.filter.pscyst <- prune_taxa(filter.pscyst.rab, pscyst)


filter.pscyst.mid.rab <- phyloseq::genefilter_sample(pscystmid.rab, filterfun_sample(function(x) x / sum(x) > .01))
forancom.filter.pscystmid <- prune_taxa(filter.pscyst.mid.rab, pscystmid)

filter.pscyst.fall.rab <- phyloseq::genefilter_sample(pscystfall.rab, filterfun_sample(function(x) x / sum(x) > .01))
forancom.filter.pscystfall <- prune_taxa(filter.pscyst.fall.rab, pscystfall)


psbulk.rab <- transform_sample_counts(psbulk, function(x) x/sum(x))
filter.psbulk.rab <- phyloseq::genefilter_sample(psbulk.rab, filterfun_sample(function(x) x / sum(x) > .02))
#pull from the non relative abundance bulk file
forancom.filter.psbulk <- prune_taxa(filter.psbulk.rab, psbulk)


filter.psbulk.mid.rab <- phyloseq::genefilter_sample(psbulkmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.psbulkmid <- prune_taxa(filter.psbulk.mid.rab, psbulkmid)


filter.psbulk.fall.rab <- phyloseq::genefilter_sample(psbulkfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.psbulkfall <-prune_taxa(filter.psbulk.fall.rab, psbulkfall)



psrhiz.rab <- transform_sample_counts(psrhiz, function(x) x/sum(x))
filter.psrhiz.rab <- phyloseq::genefilter_sample(psrhiz.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.psrhiz <- prune_taxa(filter.psrhiz.rab, psrhiz)


psroot.rab <- transform_sample_counts(psroot, function(x) x/sum(x))
filter.psroot.rab <- phyloseq::genefilter_sample(psroot.rab, filterfun_sample(function(x) x / sum(x) > .02))
forancom.filter.psroot <- prune_taxa(filter.psroot.rab, psroot)


#Ps5 has abundance cut offs for each sample! 
ps5 <- merge_phyloseq(forancom.filter.pscyst, forancom.filter.psbulk, forancom.filter.psrhiz, forancom.filter.psroot)
#check to see if it works
testancombc.ps5 <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                            group = "Sample_Type", n_cl = 8, global = TRUE)

#Bulk as comparator
sample_data(ps5)$Sample_Type <- as.factor(sample_data(ps5)$Sample_Type)
ancombc.bulk <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.bulk = ancombc.bulk$res
#res.bulk_global = ancombc.bulk$res_global


#Cyst 
sample_data(ps5)$Sample_Type <- relevel(sample_data(ps5)$Sample_Type, "Cyst")
ancombc.cyst <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.cyst = ancombc.cyst$res
#res.cyst_global = ancombc.cyst$res_global


#Rhizosphere
sample_data(ps5)$Sample_Type <- relevel(sample_data(ps5)$Sample_Type, "Rhizosphere")
ancombc.rhiz <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.rhiz = ancombc.rhiz$res
#res.rhiz_global = ancombc.rhiz$res_global


#Root
sample_data(ps5)$Sample_Type <- relevel(sample_data(ps5)$Sample_Type, "Roots")
ancombc.root <- ancombc2(ps5 ,assay_name = "counts", tax_level = "Genus", fix_formula = "Sample_Type",
                         group = "Sample_Type", n_cl = 8, global = TRUE)
res.root = ancombc.root$res
#res.root_global = ancombc.root$res_global

write.csv(res.bulk, "res.ITS.bulk.csv")
write.csv(res.cyst, "res.ITS.cyst.csv")
write.csv(res.rhiz, "res.ITS.rhiz.csv")
write.csv(res.root, "res.ITS.root.csv")


#Using the res.**.csv files.  
#Removed anything that was false across all treatments and saved as a new csv. 
#made 3 csv files that have the 4 columns. Genera, lfc T or F values, lfc number values, and the X vs Y
#I can merge these files now to make a dataframe and try a heatmap 

cvb <- read.csv("cystvsbulk.ITS.csv")
cvrh <- read.csv("cystvrhiz.ITS.csv")
cvr <- read.csv("cystvsroot.ITS.csv")

cystvsall.ancom.df <- rbind(cvb, cvrh, cvr)

ggplot(cystvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom ITS All Samples- Cyst.pdf"

#bulk
bvc <- read.csv("bulkvcyst.ITS.csv")
bvrh <- read.csv("bulkvrhiz.ITS.csv")
bvr <- read.csv("bulkvroot.ITS.csv")

bulkvsall.ancom.df <- rbind(bvc, bvrh, bvr)

ggplot(bulkvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom ITS All Samples- Bulk.pdf"

#rhiz
rhvb <- read.csv("rhizvbulk.ITS.csv")
rhvc <- read.csv("rhizvcyst.ITS.csv")
rhvr <- read.csv("rhizvroot.ITS.csv")

rhizvsall.ancom.df <- rbind(rhvb, rhvc, rhvr)

ggplot(rhizvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom ITS All Samples- Rhizz.pdf"

#roots
rvb <- read.csv("rootvbulk.ITS.csv")
rvc <- read.csv("rootvcyst.ITS.csv")
rvrh <- read.csv("rootvrhiz.ITS.csv")

rootvsall.ancom.df <- rbind(rvb, rvc, rvrh)

ggplot(rootvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom ITS All Samples- Root.pdf"

#all sample types & comparisons in one.
all.ancom.df <- rbind(cystvsall.ancom.df, bulkvsall.ancom.df, rhizvsall.ancom.df, rootvsall.ancom.df)

ggplot(all.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom ITS Entire Dataset- merged for all sample types.pdf"


## ANCOM-BC- Entire dataset, compare treatments  
### S3
sample_data(ps5)$Plot <- as.factor(sample_data(ps5)$Plot)

sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "S3")

ancombc.ps5.S3 <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.ps5.S3 = ancombc.ps5.S3$res


### Sa
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "Sa")

ancombc.ps5.Sa <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.ps5.Sa = ancombc.ps5.Sa$res


### Sr
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "Sr")

ancombc.ps5.Sr <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.ps5.Sr = ancombc.ps5.Sr$res


### Ss
sample_data(ps5)$Plot <- relevel(sample_data(ps5)$Plot, "Ss")

ancombc.ps5.Ss <- ancombc2(ps5, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.ps5.Ss = ancombc.ps5.Ss$res

write.csv(res.ps5.S3, "res.ITS.ps5.S3.csv")

write.csv(res.ps5.Sa, "res.ITS.ps5.Sa.csv")

write.csv(res.ps5.Sr, "res.ITS.ps5.Sr.csv")

write.csv(res.ps5.Ss, "res.ITS.ps5.Ss.csv")

#Using the res.ps5.S#.csv- I removed anything that was false across all treatments and saved as filt.res.ITS.ps5.S#.csv
#made 3 csv files that have the 4 columns. Genera, lfc T or F values, lfc number values, and the treatment vs treatment
#I can merge these files now to make a dataframe and try a heatmap 

### S3 results
ps5.s3vsa <- read.csv("S3vsSA.ITS.csv")

ps5.s3vsr <- read.csv("S3vsSR.ITS.csv")

ps5.s3vss <- read.csv("S3vsSS.ITS.csv")

ps5.s3vsall.ancom.df <- rbind(ps5.s3vsa, ps5.s3vsr, ps5.s3vss)

ggplot(ps5.s3vsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS All samples- S3.pdf"

### SA results 
ps5.savs3 <- read.csv("SAvsS3.ITS.csv")

ps5.savsr <- read.csv("SAvsSR.ITS.csv")

ps5.savss <- read.csv("SAvsSS.ITS.csv")

ps5.savsall.ancom.df <- rbind(ps5.savs3, ps5.savsr, ps5.savss)

ggplot(ps5.savsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS All samples- SA.pdf"

### SR results 
ps5.srvs3 <- read.csv("SRvsS3.ITS.csv")

ps5.srvsa <- read.csv("SRvsSA.ITS.csv")

ps5.srvss <- read.csv("SRvsSS.ITS.csv")

ps5.srvsall.ancom.df <- rbind(ps5.srvs3, ps5.srvsa, ps5.srvss)

ggplot(ps5.srvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS All samples- SR.pdf"

### SS results 
ps5.ssvs3 <- read.csv("SSvsS3.ITS.csv")

ps5.ssvsa <- read.csv("SSvsSA.ITS.csv")

ps5.ssvsr <- read.csv("SSvsSR.ITS.csv")

ps5.ssvsall.ancom.df <- rbind(ps5.ssvs3, ps5.ssvsa, ps5.ssvsr)

ggplot(ps5.ssvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS All samples- SS.pdf"

### merge all treatment comparisons in one.
ps5.alltreatment.ancom.df <- rbind(ps5.s3vsall.ancom.df, ps5.savsall.ancom.df, ps5.srvsall.ancom.df, ps5.ssvsall.ancom.df)

write.csv(ps5.alltreatment.ancom.df, "ps5.alltreatment.ancom.csv")

#removed all the false lfc values

ps5.alltreatment.ancom.df <- read.csv("ps5.alltreatment.ancom.csv")

ggplot(ps5.alltreatment.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS-  Entire dataset-merged to compare all treatments.pdf"
 

## ANCOM-BC- Just cysts- compare the treatments- NOTHING WAS SIGNIFICANT  
#Just cysts- compare treatments

#using the filtered phyloseq object of cysts

forancom.filter.pscyst

### S3
sample_data(forancom.filter.pscyst)$Plot <- as.factor(sample_data(forancom.filter.pscyst)$Plot)

sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "S3")

ancombc.cystS3 <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.cystS3 = ancombc.cystS3$res

write.csv(res.cystS3, "res.cystS3.csv")

### Sa
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "Sa")

ancombc.cystSA <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.cystSA = ancombc.cystSA$res

write.csv(res.cystSA, "res.cystSA.csv")

### Sr
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "Sr")

ancombc.cystSR <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)
res.cystSR = ancombc.cystSR$res

write.csv(res.cystSR, "res.cystSR.csv")

### Ss
sample_data(forancom.filter.pscyst)$Plot <- relevel(sample_data(forancom.filter.pscyst)$Plot, "Ss")

ancombc.cystSS <- ancombc2(forancom.filter.pscyst ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.cystSS = ancombc.cystSS$res

write.csv(res.cystSS, "res.cystSS.csv")

#Nothing was significant! 



## ANCOM-BC- Just mid season cysts- compare the treatments  
#using the filtered phyloseq object of cysts

forancom.filter.pscystmid 

### S3
sample_data(forancom.filter.pscystmid)$Plot <- as.factor(sample_data(forancom.filter.pscystmid)$Plot)

sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "S3")

ancombc.cystmidS3 <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.cystmidS3 = ancombc.cystmidS3$res

#write.csv(res.cystmidS3, "res.cystmidS3.csv")


### Sa
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "Sa")

ancombc.cystmidSA <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.cystmidSA = ancombc.cystmidSA$res

#write.csv(res.cystmidSA, "res.cystmidSA.csv")

### Sr
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "Sr")

ancombc.cystmidSR <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.cystmidSR = ancombc.cystmidSR$res

#write.csv(res.cystmidSR, "res.cystmidSR.csv")

### Ss
sample_data(forancom.filter.pscystmid)$Plot <- relevel(sample_data(forancom.filter.pscystmid)$Plot, "Ss")

ancombc.cystmidSS <- ancombc2(forancom.filter.pscystmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.cystmidSS = ancombc.cystmidSS$res

#write.csv(res.cystmidSS, "res.cystmidSS.csv")

### View Cyst mid ressults 
cystmid.savss <- read.csv("cystmidSAvscystmidSS.ITS.csv")

cystmid.ssvsa <- read.csv("cystmidSSvscystmidSA.ITS.csv")

ggplot(cystmid.savss, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as Ancom ITS Cyst Mid SA vs all.pdf"

ggplot(cystmid.ssvsa, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as Ancom ITS Cyst Mid SS vs all.pdf"

ancom.cystmid.df <- rbind(cystmid.savss, cystmid.ssvsa)

ggplot(ancom.cystmid.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as Ancom ITS Cyst Mid- merged treatments.pdf"


## ANCOM-BC- Just fall season cysts- compare the treatments- NOTHING WAS SIGNIFICANT 
#using the filtered phyloseq object of cysts

forancom.filter.pscystfall 

### S3
sample_data(forancom.filter.pscystfall)$Plot <- as.factor(sample_data(forancom.filter.pscystfall)$Plot)

sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "S3")

ancombc.cystfallS3 <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)


res.cystfallS3 = ancombc.cystfallS3$res

#write.csv(res.cystfallS3, "res.cystfallS3.csv")

### Sa
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "Sa")

ancombc.cystfallSA <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.cystfallSA = ancombc.cystfallSA$res

#write.csv(res.cystfallSA, "res.cystfallSA.csv")

### Sr
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "Sr")

ancombc.cystfallSR <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.cystfallSR = ancombc.cystfallSR$res

#write.csv(res.cystfallSR, "res.cystfallSR.csv")

### Ss
sample_data(forancom.filter.pscystfall)$Plot <- relevel(sample_data(forancom.filter.pscystfall)$Plot, "Ss")

ancombc.cystfallSS <- ancombc2(forancom.filter.pscystfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.cystfallSS = ancombc.cystfallSS$res

#write.csv(res.cystfallSS, "res.cystfallSS.csv")



## ANCOM-BC - Just cysts- compare season 
forancom.filter.pscyst

ancombc.cyst.season <- ancombc2(forancom.filter.pscyst, assay_name = "counts", tax_level = "Genus", fix_formula = "Season",
                                group = "Season", n_cl = 8, global = TRUE)

res.cyst.season = ancombc.cyst.season$res

#write.csv(res.cystseason, "res.cystseason.csv")

ancom.cyst.season <- read.csv("cyst.season.ITS.csv")

ggplot(ancom.cyst.season, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Cyst-season.pdf"



## ANCOM-BC- Just bulk- compare the treatments  
forancom.filter.psbulk

### S3 
sample_data(forancom.filter.psbulk)$Plot <- as.factor(sample_data(forancom.filter.psbulk)$Plot)

sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "S3")

ancombc.bulkS3 <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.bulkS3 = ancombc.bulkS3$res

write.csv(res.bulkS3, "res.bulkS3.csv")

### Sa
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "Sa")

ancombc.bulkSA <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.bulkSA = ancombc.bulkSA$res

write.csv(res.bulkSA, "res.bulkSA.csv")

### Sr
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "Sr")

ancombc.bulkSR <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.bulkSR = ancombc.bulkSR$res

write.csv(res.bulkSR, "res.bulkSR.csv")

### Ss
sample_data(forancom.filter.psbulk)$Plot <- relevel(sample_data(forancom.filter.psbulk)$Plot, "Ss")

ancombc.bulkSS <- ancombc2(forancom.filter.psbulk ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.bulkSS = ancombc.bulkSS$res

write.csv(res.bulkSS, "res.bulkSS.csv")

### view bulk soil results 
bulk.savs3 <- read.csv("bulkSAvbulkS3.ITS.csv")

bulk.savsr <- read.csv("bulkSAvbulkSR.ITS.csv")

bulk.savss <- read.csv("bulkSAvbulkSS.ITS.csv")

bulk.srvs3 <- read.csv("bulkSRvbulkS3.ITS.csv")

bulk.srvsa <- read.csv("bulkSRvbulkSA.ITS.csv")

bulk.srvss <- read.csv("bulkSRvbulkSS.ITS.csv")

bulk.ssvs3 <- read.csv("bulkSSvbulkS3.ITS.csv")

bulk.ssvsa <- read.csv("bulkSSvbulkSA.ITS.csv")

bulk.ssvsr <- read.csv("bulkSSvbulkSR.ITS.csv")


bulk.savsall.ancom.df <- rbind(bulk.savs3, bulk.savsr, bulk.savss)

ggplot(bulk.savsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk SA vs All.pdf"


bulk.srvsall.ancom.df <- rbind(bulk.srvs3, bulk.srvsa, bulk.srvss)

ggplot(bulk.srvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk SR vs All.pdf"

bulk.ssvsall.ancom.df <- rbind(bulk.ssvs3, bulk.ssvsa, bulk.ssvsr)

ggplot(bulk.ssvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk SS vs All.pdf"

### merge bulk soil samples together 
ancom.bulk.df <- rbind(bulk.savsall.ancom.df, bulk.srvsall.ancom.df, bulk.ssvsall.ancom.df)

ggplot(ancom.bulk.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk merged treatments.pdf"



## ANCOM-BC- Just mid bulk- compare the treatments 
forancom.filter.psbulkmid
### S3
sample_data(forancom.filter.psbulkmid)$Plot <- as.factor(sample_data(forancom.filter.psbulkmid)$Plot)

sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "S3")

ancombc.bulkmidS3 <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)


res.bulkmidS3 = ancombc.bulkmidS3$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidS3, "res.bulkmidS3.csv")

### Sa
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "Sa")

ancombc.bulkmidSA <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.bulkmidSA = ancombc.bulkmidSA$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidSA, "res.bulkmidSA.csv")

### Sr
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "Sr")

ancombc.bulkmidSR <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)

res.bulkmidSR = ancombc.bulkmidSR$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidSR, "res.bulkmidSR.csv")

### Ss
sample_data(forancom.filter.psbulkmid)$Plot <- relevel(sample_data(forancom.filter.psbulkmid)$Plot, "Ss")

ancombc.bulkmidSS <- ancombc2(forancom.filter.psbulkmid ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                              group = "Plot", n_cl = 8, global = TRUE)


res.bulkmidSS = ancombc.bulkmidSS$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkmidSS, "res.bulkmidSS.csv")


### view bulk mid results 
bulkmid.s3vsr <- read.csv("bulkmidS3vbulkmidSR.ITS.csv")

bulkmid.s3vss <- read.csv("bulkmidS3vbulkmidSS.ITS.csv")

bulkmid.savsr <- read.csv("bulkmidSAvbulkmidSR.ITS.csv")

bulkmid.savss <- read.csv("bulkmidSAvbulkmidSS.ITS.csv")

bulkmid.srvs3 <- read.csv("bulkmidSRvbulkmidS3.ITS.csv")

bulkmid.ssvs3 <- read.csv("bulkmidSSvbulkmidS3.ITS.csv")

bulkmid.ssvsa <- read.csv("bulkmidSSvbulkmidSA.ITS.csv")


ancom.bulkmid.s3vall.df <- rbind(bulkmid.s3vsr, bulkmid.s3vss)

ggplot(ancom.bulkmid.s3vall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")
#Saved as "Ancom ITS Bulk Mid S3 vs all.pdf"

ancom.bulkmid.savall.df <- rbind(bulkmid.savsr, bulkmid.savss)

ggplot(ancom.bulkmid.savall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Mid SA vs all.pdf"


ggplot(bulkmid.srvs3, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Mid SR vs all.pdf"

ancom.bulkmid.ssvall.df <- rbind(bulkmid.ssvs3, bulkmid.ssvsa)

ggplot(ancom.bulkmid.ssvall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Mid SS vs all.pdf"

### merge bulk mid results together 
ancom.bulkmid.df <- rbind(ancom.bulkmid.s3vall.df, ancom.bulkmid.savall.df, bulkmid.srvs3, ancom.bulkmid.ssvall.df)

ggplot(ancom.bulkmid.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Mid- merged treatments.pdf"

 

## ANCOM-BC- Just fall bulk- compare the treatments  
forancom.filter.psbulkfall 
### S3
sample_data(forancom.filter.psbulkfall)$Plot <- as.factor(sample_data(forancom.filter.psbulkfall)$Plot)

sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "S3")

ancombc.bulkfallS3 <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.bulkfallS3 = ancombc.bulkfallS3$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkfallS3, "res.bulkfallS3.csv")

### Sa
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "Sa")

ancombc.bulkfallSA <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.bulkfallSA = ancombc.bulkfallSA$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkfallSA, "res.bulkfallSA.csv")

### Sr
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "Sr")

ancombc.bulkfallSR <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)

res.bulkfallSR = ancombc.bulkfallSR$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkfallSR, "res.bulkfallSR.csv")

### Ss
sample_data(forancom.filter.psbulkfall)$Plot <- relevel(sample_data(forancom.filter.psbulkfall)$Plot, "Ss")

ancombc.bulkfallSS <- ancombc2(forancom.filter.psbulkfall ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                               group = "Plot", n_cl = 8, global = TRUE)


res.bulkfallSS = ancombc.bulkfallSS$res

#Copied from dataframe into an excel sheet to save time. write.csv(res.bulkfallSS, "res.bulkfallSS.csv")

### View bulk fall results 
bulkfall.s3vsr <-read.csv("bulkfallS3vbulkfallSR.ITS.csv")

bulkfall.s3vss <-read.csv("bulkfallS3vbulkfallSS.ITS.csv")

bulkfall.savsr <-read.csv("bulkfallSAvbulkfallSR.ITS.csv")

bulkfall.savss <-read.csv("bulkfallSAvbulkfallSS.ITS.csv")

bulkfall.srvs3 <- read.csv("bulkfallSRvbulkfallS3.ITS.csv")

bulkfall.srvsa <-read.csv("bulkfallSRvbulkfallSA.ITS.csv")

bulkfall.ssvs3 <-read.csv("bulkfallSSvbulkfallS3.ITS.csv")

bulkfall.ssvsa <-read.csv("bulkfallSSvbulkfallSA.ITS.csv")


ancom.bulkfall.s3vall.df <- rbind(bulkfall.s3vsr, bulkfall.s3vss)

ggplot(ancom.bulkfall.s3vall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Fall S3 vs all.pdf"

ancom.bulkfall.savall.df <- rbind(bulkfall.savsr, bulkfall.savss)

ggplot(ancom.bulkfall.savall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Fall SA vs all.pdf"

ancom.bulkfall.srvall.df <- rbind(bulkfall.srvs3, bulkfall.srvsa)

ggplot(ancom.bulkfall.srvall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Fall SR vs all.pdf"

ancom.bulkfall.ssvall.df <- rbind(bulkfall.ssvs3, bulkfall.ssvsa)

ggplot(ancom.bulkfall.ssvall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Fall SS vs all.pdf"

### merge all bulk fall results together
ancom.allbulkfall.df <- rbind(ancom.bulkfall.s3vall.df, ancom.bulkfall.savall.df, ancom.bulkfall.srvall.df, ancom.bulkfall.ssvall.df)

ggplot(ancom.allbulkfall.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk Fall- all treatments.pdf"


## ANCOM-BC - Just bulk- compare season 
forancom.filter.psbulk

ancombc.bulk.season <- ancombc2(forancom.filter.psbulk, assay_name = "counts", tax_level = "Genus", fix_formula = "Season",
                                group = "Season", n_cl = 8, global = TRUE)

res.bulk.season = ancombc.bulk.season$res

ancom.bulk.season <- read.csv("bulk.season.ITS.csv")

ggplot(ancom.bulk.season, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Bulk-season.pdf"

## ANCOM-BC- Just rhiz- compare the treatments  
forancom.filter.psrhiz
### S3
sample_data(forancom.filter.psrhiz)$Plot <- as.factor(sample_data(forancom.filter.psrhiz)$Plot)

sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "S3")

ancombc.rhizS3 <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rhizS3 = ancombc.rhizS3$res

write.csv(res.rhizS3, "res.rhizS3.csv")

### Sa
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "Sa")

ancombc.rhizSA <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rhizSA = ancombc.rhizSA$res

write.csv(res.rhizSA, "res.rhizSA.csv")

### Sr
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "Sr")

ancombc.rhizSR <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rhizSR = ancombc.rhizSR$res

write.csv(res.rhizSR, "res.rhizSR.csv")

### Ss
sample_data(forancom.filter.psrhiz)$Plot <- relevel(sample_data(forancom.filter.psrhiz)$Plot, "Ss")

ancombc.rhizSS <- ancombc2(forancom.filter.psrhiz ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rhizSS = ancombc.rhizSS$res

write.csv(res.rhizSS, "res.rhizSS.csv")

### View Rhizosphere ancom results 
rhiz.s3vsa <- read.csv("rhizS3vsrhizSA.ITS.csv")

rhiz.s3vsr <- read.csv("rhizS3vsrhizSR.ITS.csv")

rhiz.s3vss <- read.csv("rhizS3vsrhizSS.ITS.csv")

rhiz.savs3 <- read.csv("rhizSAvsrhizS3.ITS.csv")

rhiz.savsr <- read.csv("rhizSAvsrhizSR.ITS.csv")

rhiz.savss <- read.csv("rhizSAvsrhizSS.ITS.csv")

rhiz.srvs3 <- read.csv("rhizSRvsrhizS3.ITS.csv")

rhiz.srvsa <- read.csv("rhizSRvsrhizSA.ITS.csv")

rhiz.srvss <- read.csv("rhizSRvsrhizSS.ITS.csv")


rhiz.s3vsall.ancom.df <- rbind(rhiz.s3vsa, rhiz.s3vsr, rhiz.s3vss)

ggplot(rhiz.s3vsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Rhiz S3 vs All.pdf"


rhiz.savsall.ancom.df <- rbind(rhiz.savs3, rhiz.savsr, rhiz.savss)

ggplot(rhiz.savsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Rhiz SA vs All.pdf"

rhiz.srvsall.ancom.df <- rbind(rhiz.srvs3, rhiz.srvsa, rhiz.srvss)

ggplot(rhiz.srvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Rhiz SR vs All.pdf"

### merge all rhizosphere results together
ancom.rhiz.df <- rbind(rhiz.s3vsall.ancom.df, rhiz.savsall.ancom.df, rhiz.srvsall.ancom.df)

write.csv(ancom.rhiz.df, "ancom.rhiz.csv")

#remove false lfc values

ancom.rhiz.df <- read.csv("ancom.rhiz.csv")

ggplot(ancom.rhiz.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Rhiz- merged treatments.pdf"


## ANCOM-BC- Just roots- compare the treatments  
forancom.filter.psroot
### S3 
sample_data(forancom.filter.psroot)$Plot <- as.factor(sample_data(forancom.filter.psroot)$Plot)

sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "S3")

ancombc.rootS3 <- ancombc2(forancom.filter.psroot ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rootS3 = ancombc.rootS3$res

write.csv(res.rootS3, "res.rootsS3.csv")

### Sa
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "Sa")

ancombc.rootSA <- ancombc2(forancom.filter.psroot ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rootSA = ancombc.rootSA$res

write.csv(res.rootSA, "res.rootsSA.csv")

### Sr
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "Sr")

ancombc.rootSR <- ancombc2(forancom.filter.psroot, assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rootSR = ancombc.rootSR$res

write.csv(res.rootSR, "res.rootSR.csv")

### Ss
sample_data(forancom.filter.psroot)$Plot <- relevel(sample_data(forancom.filter.psroot)$Plot, "Ss")

ancombc.rootSS <- ancombc2(forancom.filter.psroot ,assay_name = "counts", tax_level = "Genus", fix_formula = "Plot",
                           group = "Plot", n_cl = 8, global = TRUE)

res.rootSS = ancombc.rootSS$res

write.csv(res.rootSS, "res.rootsSS.csv")

## Root ancom results 
root.s3vsa <- read.csv("rootS3vsrootSA.ITS.csv")

root.s3vsr <- read.csv("rootS3vsrootSR.ITS.csv")

root.s3vss <- read.csv("rootS3vsrootSS.ITS.csv")

root.savs3 <- read.csv("rootSAvsrootS3.ITS.csv")

root.savsr <- read.csv("rootSAvsrootSR.ITS.csv")

root.savss <- read.csv("rootSAvsrootSS.ITS.csv")

root.srvs3 <- read.csv("rootSRvsrootS3.ITS.csv")

root.srvsa <- read.csv("rootSRvsrootSA.ITS.csv")

root.srvss <- read.csv("rootSRvsrootSS.ITS.csv")

root.ssvs3 <- read.csv("rootSSvsrootS3.ITS.csv")

root.ssvsa <- read.csv("rootSSvsrootSA.ITS.csv")

root.ssvsr <- read.csv("rootSSvsrootSR.ITS.csv")


root.s3vsall.ancom.df <- rbind(root.s3vsa, root.s3vsr, root.savss)

ggplot(root.s3vsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Root S3 vs All.pdf"

root.savsall.ancom.df <- rbind(root.savs3, root.savsr, root.savss)

ggplot(root.savsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Root SA vs All.pdf"


root.srvsall.ancom.df <- rbind(root.srvs3, root.srvsa, root.srvss)

ggplot(root.srvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Root SR vs All.pdf"

root.ssvsall.ancom.df <- rbind(root.ssvs3, root.ssvsa, root.ssvsr)

ggplot(root.ssvsall.ancom.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")


#Saved as "Ancom ITS Root SS vs All.pdf"

### merge all root results together
ancom.root.df <- rbind(root.s3vsall.ancom.df, root.savsall.ancom.df, root.srvsall.ancom.df, root.ssvsall.ancom.df)

write.csv(ancom.root.df, "ancom.root.csv")

#remove false lfc values

ancom.root.df <- read.csv("ancom.root.csv")

ggplot(ancom.root.df, aes(x= lfc_TF,y= Genus, fill = lfc_value)) + geom_tile(color = "black", 
) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey"), axis.text.y.left = element_text(size = 12),
          strip.text.x.top = element_text(size = 10)) + facet_grid(~Sample_vs, scales = "free", space = "free"
          ) + scale_fill_gradient2(low = "red", mid = "white", high = "blue")

#Saved as "Ancom ITS Root- merged treatments.pdf"
