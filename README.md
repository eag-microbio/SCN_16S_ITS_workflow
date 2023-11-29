# SCN 16S
---
title: "LTR16S_S3ARS_Markdown"
author: "Emily Green"
date: "2023-11-02"
output: html_document
---
## Import data & assemble

### Run DADA2

```{r dada2, echo=FALSE, include=FALSE, cache=TRUE, results = "hide"}
library(dada2)
path <- "/workdir/eag252/LTR_16S_S3ARS/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
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

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
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

taxa <- assignTaxonomy(seqtab.nochim, "/workdir/eag252/LTR_16S_S3ARS/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

#write.csv(track, "all16sfiltcounts.csv")
#Output the taxa, taxa counts, and taxa asv sequences
#write.csv(taxa.print, "taxaLTR16Sfilt.csv")
#write.csv(seqtab.nochim, "ASVLTR16Sfilt.csv")
```

### Load packages
```{r loadpackages, echo=FALSE, include=FALSE}
library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)
library(viridisLite)
library(viridis)
library(ape)
library(cowplot)
theme_set(theme_bw())
```

### Import metadata
```{r importmetadata, echo=FALSE}
meta <- read.csv("/workdir/eag252/LTR_16S_S3ARS/LTR16s_S3ARSMetadata.csv",header = TRUE, row.names = 1)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(meta), tax_table(taxa))
#sample_data(ps)
#ps
```

## Quality Control

### Filter taxa
```{r filtertaxa, echo=FALSE, include=FALSE, results = "hide"}
sort(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))

ps1 <- subset_taxa(ps,Kingdom == "Bacteria")
get_taxa_unique(ps1, taxonomic.rank = "Kingdom")
ps1

sort(get_taxa_unique(ps1, taxonomic.rank = "Phylum"))
sort(get_taxa_unique(ps1, taxonomic.rank = "Family"))
sort(get_taxa_unique(ps1, taxonomic.rank = "Order"))

#Remove mitrochondria & Chloroplast
ps2 <- subset_taxa(ps1, !Family == "Mitochondria")
sort(get_taxa_unique(ps2, taxonomic.rank = "Family"))
ps2 <- subset_taxa(ps2, !Order == "Chloroplast")
ps2
```

### Quick check alpha diversity
```{r quickalpha, echo=FALSE}
plot_richness(ps2, x="SampleType", measures=c("Shannon", "Simpson"), color="Sample_or_Control")
```

### Decontam
```{r decontam, echo=FALSE, message=FALSE}
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

#sample_sums(ps2)
#View(sample_sums(ps2))
```

```{r removesampleslessthan100, echo=FALSE}
ps3 <- prune_samples(sample_sums(ps2) > 100, ps2)
#View(sample_sums(ps3))
```

### Quick Unifrac to check all samples & controls
```{r echo=FALSE}
ps3t <- ps3

ps3t.tree <- rtree(ntaxa(ps3t), rooted = TRUE, tip.label = taxa_names(ps3t))
#plot(ps3t.tree)
ps3t <- merge_phyloseq(ps3t, ps3t.tree)
```

#### Unweighted
```{r echo=FALSE}
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=F)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "Sample_or_Control", shape = "SampleType") + theme(aspect.ratio=1)
```

```{r}
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~SampleType*Sample_or_Control, ps3t.tree.data)
```

#### Weighted
```{r echo=FALSE}
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=T)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "Sample_or_Control", shape = "SampleType") + theme(aspect.ratio=1)
```

```{r}
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~SampleType*Sample_or_Control, ps3t.tree.data)
```

### Remove all negative controls and mock communities
```{r removecontrols, echo=FALSE}
ps4 <- subset_samples(ps3, Sample_or_Control == "Sample")
#View(sample_data(ps4))
```

### Quick view of beta-diversity of soil, cysts and roots.

```{r echo=FALSE}
ps4t <-ps4 

#Unifrac
ps4t.tree <- rtree(ntaxa(ps4t), rooted = TRUE, tip.label = taxa_names(ps4t))
#plot(ps4t.tree)
ps4t <- merge_phyloseq(ps4t, ps4t.tree)
```

#### Unweighted
```{r echo=FALSE}
ps4tunifrac_dist <- phyloseq::distance(ps4t, method="unifrac", weighted=F)
ordination = ordinate(ps4t, method="PCoA", distance=ps4tunifrac_dist)
plot_ordination(ps4t, ordination, color = "SampleType", shape = "Season") + theme(aspect.ratio=1)
```

```{r}
ps4t.tree.data <- data.frame(sample_data(ps4t))
adonis2(ps4tunifrac_dist ~SampleType*Season, ps4t.tree.data)
```

#### Weighted
```{r echo=FALSE}
ps4tunifrac_dist <- phyloseq::distance(ps4t, method="unifrac", weighted=T)
ordination = ordinate(ps4t, method="PCoA", distance=ps4tunifrac_dist)
plot_ordination(ps4t, ordination, color = "SampleType", shape = "Season") + theme(aspect.ratio=1)
```

```{r}
ps4t.tree.data <- data.frame(sample_data(ps4t))
adonis2(ps4tunifrac_dist ~SampleType*Season, ps4t.tree.data)
```

## Split samples based on Sample Type into separate phyloseq objects
```{r splitsamplestophyloseq}
pssoil <- subset_samples(ps4, SampleType == "Soil")
pscyst <- subset_samples(ps4, SampleType == "Cyst")
psroot <- subset_samples(ps4, SampleType == "Roots")
```

### Confirm that each phyloseq object contains the correct samples

```{r}
#sample_data(pssoil)
#sample_data(pscyst)
#sample_data(psroot)
```

### Remove taxa with 0 reads in each phyloseq object
```{r remove0readtaxa}
pssoil <- prune_taxa(taxa_sums(pssoil) > 0, pssoil)
pscyst <- prune_taxa(taxa_sums(pscyst) > 0, pscyst)
psroot <- prune_taxa(taxa_sums(psroot) > 0, psroot)
```


## Alpha Diversity & Anova statistics

### Soil Alpha Diversity
```{r soilalpha, echo=FALSE, message = FALSE}
soilalpha <- estimate_richness(pssoil, measures = c("Shannon", "Simpson"))
soilalpha$Plot <- sample_data(pssoil)$Plot
soilalpha$SoilType <- sample_data(pssoil)$SoilType
soilalpha$Season <- sample_data(pssoil)$Season
```

#### Soil Shannon Alpha diversity
```{r echo=FALSE}
ggplot(soilalpha, aes(x = Plot, y = Shannon, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
```

```{r}
soilalpha.soiltype <- aov(Shannon ~SoilType, soilalpha)
anova(soilalpha.soiltype)
soilalpha.season <- aov(Shannon ~Season, soilalpha)
anova(soilalpha.season)
soilalpha.plot <- aov(Shannon ~Plot, soilalpha)
anova(soilalpha.plot)
```

#### Soil Simpson Alpha Diversity
```{r echo=FALSE}
ggplot(soilalpha, aes(x = Plot, y = Simpson, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
```

```{r}
soilalpha.soiltype <- aov(Simpson ~SoilType, soilalpha)
anova(soilalpha.soiltype)
soilalpha.season <- aov(Simpson ~Season, soilalpha)
anova(soilalpha.season)
soilalpha.plot <- aov(Simpson ~Plot, soilalpha)
anova(soilalpha.plot)
```

### Cysts Alpha Diversity
```{r cystalpha, message = FALSE, echo=FALSE}
cystalpha <- estimate_richness(pscyst,measures = c("Shannon", "Simpson"))
cystalpha$Plot <- sample_data(pscyst)$Plot
cystalpha$Season <- sample_data(pscyst)$Season
```

#### Cyst Shannon Alpha diversity
```{r echo=FALSE}
ggplot(cystalpha, aes(x= Plot, y = Shannon, color = Season)) + geom_boxplot()
```

```{r}
cystalpha.plot <- aov(Shannon ~Plot, cystalpha)
anova(cystalpha.plot)
cystalpha.season <- aov(Shannon ~Season, cystalpha)
anova(cystalpha.season)
```

#### Cyst Simpson Alpha Diversity
```{r echo=FALSE}
ggplot(cystalpha, aes(x= Plot, y = Simpson, color = Season)) + geom_boxplot()
```

```{r}
cystalpha.plot <- aov(Simpson ~Plot, cystalpha)
anova(cystalpha.plot)
cystalpha.season <- aov(Simpson ~Season, cystalpha)
anova(cystalpha.season)
```

### Roots Alpha Diversity
```{r rootalpha, message = FALSE, echo=FALSE}
rootalpha <- estimate_richness(psroot, measures = c("Shannon", "Simpson"))
rootalpha$Plot <- sample_data(psroot)$Plot
```

#### Roots Shannon Alpha diversity
```{r echo=FALSE}
ggplot(rootalpha, aes(x= Plot, y = Shannon)) + geom_boxplot()
```

```{r}
rootalpha.plot <- aov(Shannon ~Plot, rootalpha)
anova(rootalpha.plot)
```

#### Roots Simpson Alpha Diversity
```{r echo=FALSE}
ggplot(rootalpha, aes(x= Plot, y = Simpson)) + geom_boxplot()
```

```{r}
rootalpha.plot <- aov(Simpson ~Plot, rootalpha)
anova(rootalpha.plot)
```

## Unifrac phyloseq objects & Split sample types into seasons
```{r}
pssoil.t <- pssoil
pscyst.t <- pscyst
psroot.t <- psroot

psrhiz <- subset_samples(pssoil, SoilType == "Rhizosphere")
psrhiz.t <- psrhiz
psbulk <- subset_samples(pssoil, SoilType == "Bulk soil")
psbulk.t <- psbulk

psbulkmid <- subset_samples(psbulk, Season == "Mid")
psbulkmid.t <- psbulkmid
psbulkfall <- subset_samples(psbulk, Season == "Fall")
psbulkfall.t <- psbulkfall

pscystmid <- subset_samples(pscyst, Season == "Mid")
pscystmid.t <- pscystmid
pscystfall <- subset_samples(pscyst, Season == "Fall")
pscystfall.t <- pscystfall
```

## Beta Diversity & PERMANOVA

### Soil Beta Diversity
```{r soilbeta, echo=FALSE}
pssoil.t.tree <- rtree(ntaxa(pssoil.t), rooted = TRUE, tip.label = taxa_names(pssoil.t))
#plot(pssoil.t.tree)
pssoil.t <- merge_phyloseq(pssoil.t, pssoil.t.tree)
```

#### Soil Unweighted
```{r soilunweight, echo=FALSE}
pssoiltunifrac_dist <- phyloseq::distance(pssoil.t, method="unifrac", weighted=F)
ordination = ordinate(pssoil.t, method="PCoA", distance=pssoiltunifrac_dist)
plot_ordination(pssoil.t, ordination, color = "SoilType", shape = "Season", label = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.t, ordination, color = "Plot", shape = "SoilType", label = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pssoil.t.tree.data <- data.frame(sample_data(pssoil.t))
adonis2(pssoiltunifrac_dist ~SoilType*Plot*Season, pssoil.t.tree.data)
```

#### Soil Weighted
```{r soilweight, echo=FALSE}
pssoiltunifrac_dist <- phyloseq::distance(pssoil.t, method="unifrac", weighted=T)
ordination = ordinate(pssoil.t, method="PCoA", distance=pssoiltunifrac_dist)
plot_ordination(pssoil.t, ordination, color = "SoilType", shape = "Season", label = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.t, ordination, color = "Plot", shape = "SoilType", label = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pssoil.t.tree.data <- data.frame(sample_data(pssoil.t))
adonis2(pssoiltunifrac_dist ~SoilType*Plot*Season, pssoil.t.tree.data)
```

### Rhizosphere Beta Diversity
```{r rhizbeta, echo=FALSE}
psrhiz.t.tree <- rtree(ntaxa(psrhiz.t), rooted = TRUE, tip.label = taxa_names(psrhiz.t))
#plot(psrhiz.t.tree)
psrhiz.t <- merge_phyloseq(psrhiz.t, psrhiz.t.tree)
```

#### Rhizosphere Unweighted
```{r rhizunweight, echo=FALSE}
psrhiztunifrac_dist <- phyloseq::distance(psrhiz.t, method="unifrac", weighted=F)
ordination = ordinate(psrhiz.t, method="PCoA", distance=psrhiztunifrac_dist)
plot_ordination(psrhiz.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psrhiz.t.tree.data <- data.frame(sample_data(psrhiz.t))
adonis2(psrhiztunifrac_dist ~Plot, psrhiz.t.tree.data)
```

#### Rhizosphere Weighted
```{r rhizweight, echo=FALSE}
psrhiztunifrac_dist <- phyloseq::distance(psrhiz.t, method="unifrac", weighted=T)
ordination = ordinate(psrhiz.t, method="PCoA", distance=psrhiztunifrac_dist)
plot_ordination(psrhiz.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psrhiz.t.tree.data <- data.frame(sample_data(psrhiz.t))
adonis2(psrhiztunifrac_dist ~Plot, psrhiz.t.tree.data)
```


### Bulk Soil Beta Diversity
```{r bulkbeta, echo=FALSE}
psbulk.t.tree <- rtree(ntaxa(psbulk.t), rooted = TRUE, tip.label = taxa_names(psbulk.t))
#plot(psbulk.t.tree)
psbulk.t <- merge_phyloseq(psbulk.t, psbulk.t.tree)
```

#### Bulk Unweighted
```{r bulkunweight, echo=FALSE}
psbulktunifrac_dist <- phyloseq::distance(psbulk.t, method="unifrac", weighted=F)
ordination = ordinate(psbulk.t, method="PCoA", distance=psbulktunifrac_dist)
plot_ordination(psbulk.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulk.t.tree.data <- data.frame(sample_data(psbulk.t))
adonis2(psbulktunifrac_dist ~Plot*Season, psbulk.t.tree.data)
```

#### Bulk Weighted
```{r bulkweight, echo=FALSE}
psbulktunifrac_dist <- phyloseq::distance(psbulk.t, method="unifrac", weighted=T)
ordination = ordinate(psbulk.t, method="PCoA", distance=psbulktunifrac_dist)
plot_ordination(psbulk.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulk.t.tree.data <- data.frame(sample_data(psbulk.t))
adonis2(psbulktunifrac_dist ~Plot*Season, psbulk.t.tree.data)
```

### Bulk soil Mid season Beta Diversity
```{r bulkmidbeta, echo=FALSE}
psbulkmid.t.tree <- rtree(ntaxa(psbulkmid.t), rooted = TRUE, tip.label = taxa_names(psbulkmid.t))
#plot(psbulkmid.t.tree)
psbulkmid.t <- merge_phyloseq(psbulkmid.t, psbulkmid.t.tree)
```

#### Unweighted Bulk Mid
```{r bulkmidunweight, echo=FALSE}
psbulkmidtunifrac_dist <- phyloseq::distance(psbulkmid.t, method="unifrac", weighted=F)
ordination = ordinate(psbulkmid.t, method="PCoA", distance=psbulkmidtunifrac_dist)
plot_ordination(psbulkmid.t, ordination, color = "Plot",) + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkmid.t.tree.data <- data.frame(sample_data(psbulkmid.t))
adonis2(psbulkmidtunifrac_dist ~Plot, psbulkmid.t.tree.data)
```

#### Weighted Bulk Mid
```{r bulkmidweight, echo=FALSE}
psbulkmidtunifrac_dist <- phyloseq::distance(psbulkmid.t, method="unifrac", weighted=T)
ordination = ordinate(psbulkmid.t, method="PCoA", distance=psbulkmidtunifrac_dist)
plot_ordination(psbulkmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkmid.t.tree.data <- data.frame(sample_data(psbulkmid.t))
adonis2(psbulkmidtunifrac_dist ~Plot, psbulkmid.t.tree.data)
```

### Bulk soil Fall season Beta Diversity
```{r bulkfallbeta, echo=FALSE}
psbulkfall.t.tree <- rtree(ntaxa(psbulkfall.t), rooted = TRUE, tip.label = taxa_names(psbulkfall.t))
#plot(psbulkfall.t.tree)
psbulkfall.t <- merge_phyloseq(psbulkfall.t, psbulkfall.t.tree)
```

#### Unweighted Bulk Fall
```{r bulkfallunweight, echo=FALSE}
psbulkfalltunifrac_dist <- phyloseq::distance(psbulkfall.t, method="unifrac", weighted=F)
ordination = ordinate(psbulkfall.t, method="PCoA", distance=psbulkfalltunifrac_dist)
plot_ordination(psbulkfall.t, ordination, color = "Plot",) + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkfall.t.tree.data <- data.frame(sample_data(psbulkfall.t))
adonis2(psbulkfalltunifrac_dist ~Plot, psbulkfall.t.tree.data)
```

#### Weighted Bulk Fall
```{r bulkfallweight, echo=FALSE}
psbulkfalltunifrac_dist <- phyloseq::distance(psbulkfall.t, method="unifrac", weighted=T)
ordination = ordinate(psbulkfall.t, method="PCoA", distance=psbulkfalltunifrac_dist)
plot_ordination(psbulkfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkfall.t.tree.data <- data.frame(sample_data(psbulkfall.t))
adonis2(psbulkfalltunifrac_dist ~Plot, psbulkfall.t.tree.data)
```

### Cyst Beta Diversity
```{r cystbeta, echo=FALSE}
pscyst.t.tree <- rtree(ntaxa(pscyst.t), rooted = TRUE, tip.label = taxa_names(pscyst.t))
#plot(pscyst.t.tree)
pscyst.t <- merge_phyloseq(pscyst.t, pscyst.t.tree)
```

#### Cysts Unweighted
```{r cystunweight, echo=FALSE}
pscysttunifrac_dist <- phyloseq::distance(pscyst.t, method="unifrac", weighted=F)
ordination = ordinate(pscyst.t, method="PCoA", distance=pscysttunifrac_dist)
#plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscyst.t.tree.data <- data.frame(sample_data(pscyst.t))
adonis2(pscysttunifrac_dist ~Plot*Season*Num_Cysts, pscyst.t.tree.data)
```

#### Cysts Weighted
```{r cystweight, echo=FALSE}
pscysttunifrac_dist <- phyloseq::distance(pscyst.t, method="unifrac", weighted=T)
ordination = ordinate(pscyst.t, method="PCoA", distance=pscysttunifrac_dist)
#plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscyst.t.tree.data <- data.frame(sample_data(pscyst.t))
adonis2(pscysttunifrac_dist ~Plot*Season*Num_Cysts, pscyst.t.tree.data)
```


### Cyst Mid season Beta Diversity
```{r cystmidbeta, echo=FALSE}
pscystmid.t.tree <- rtree(ntaxa(pscystmid.t), rooted = TRUE, tip.label = taxa_names(pscystmid.t))
#plot(pscystmid.t.tree)
pscystmid.t <- merge_phyloseq(pscystmid.t, pscystmid.t.tree)
```

#### Unweighted Cysts Mid
```{r cystmidunweight, echo=FALSE}
pscystmidtunifrac_dist <- phyloseq::distance(pscystmid.t, method="unifrac", weighted=F)
ordination = ordinate(pscystmid.t, method="PCoA", distance=pscystmidtunifrac_dist)
#plot_ordination(pscystmid.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystmid.t.tree.data <- data.frame(sample_data(pscystmid.t))
adonis2(pscystmidtunifrac_dist ~Plot*Num_Cysts, pscystmid.t.tree.data)
```

#### Weighted Cysts Mid
```{r cystmidweighted, echo=FALSE}
pscystmidtunifrac_dist <- phyloseq::distance(pscystmid.t, method="unifrac", weighted=T)
ordination = ordinate(pscystmid.t, method="PCoA", distance=pscystmidtunifrac_dist)
#plot_ordination(pscystmid.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystmid.t.tree.data <- data.frame(sample_data(pscystmid.t))
adonis2(pscystmidtunifrac_dist ~Plot*Num_Cysts, pscystmid.t.tree.data)
```

### Cyst Fall season Beta Diversity
```{r cystfallbeta, echo=FALSE}
pscystfall.t.tree <- rtree(ntaxa(pscystfall.t), rooted = TRUE, tip.label = taxa_names(pscystfall.t))
#plot(pscystfall.t.tree)
pscystfall.t <- merge_phyloseq(pscystfall.t, pscystfall.t.tree)
```

#### Unweighted Cysts Fall
```{r cystfallunweight, echo=FALSE}
pscystfalltunifrac_dist <- phyloseq::distance(pscystfall.t, method="unifrac", weighted=F)
ordination = ordinate(pscystfall.t, method="PCoA", distance=pscystfalltunifrac_dist)
#plot_ordination(pscystfall.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystfall.t.tree.data <- data.frame(sample_data(pscystfall.t))
adonis2(pscystfalltunifrac_dist ~Plot, pscystfall.t.tree.data)
```

#### Weighted Cysts Fall
```{r cystfallweight, echo=FALSE}
pscystfalltunifrac_dist <- phyloseq::distance(pscystfall.t, method="unifrac", weighted=T)
ordination = ordinate(pscystfall.t, method="PCoA", distance=pscystfalltunifrac_dist)
#plot_ordination(pscystfall.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystfall.t.tree.data <- data.frame(sample_data(pscystfall.t))
adonis2(pscystfalltunifrac_dist ~Plot, pscystfall.t.tree.data)
```

### Roots Beta Diversity
```{r rootbeta, echo=FALSE}
psroot.t.tree <- rtree(ntaxa(psroot.t), rooted = TRUE, tip.label = taxa_names(psroot.t))
#plot(psroot.t.tree)
psroot.t <- merge_phyloseq(psroot.t, psroot.t.tree)
```

#### Roots Unweighted
```{r rootunweight, echo=FALSE}
psroottunifrac_dist <- phyloseq::distance(psroot.t, method="unifrac", weighted=F)
ordination = ordinate(psroot.t, method="PCoA", distance=psroottunifrac_dist)
plot_ordination(psroot.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psroot.t.tree.data <- data.frame(sample_data(psroot.t))
adonis2(psroottunifrac_dist ~Plot, psroot.t.tree.data)
```

#### Roots Weighted
```{r rootweight, echo=FALSE}
psroottunifrac_dist <- phyloseq::distance(psroot.t, method="unifrac", weighted=T)
ordination = ordinate(psroot.t, method="PCoA", distance=psroottunifrac_dist)
plot_ordination(psroot.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psroot.t.tree.data <- data.frame(sample_data(psroot.t))
adonis2(psroottunifrac_dist ~Plot, psroot.t.tree.data)
```


## Bar plots highlighting the most abundant taxa
### Rhizosphere Phyla present with at least 1% abundance
```{r rhizphyla, echo=FALSE}
psrhiz.rab <- transform_sample_counts(psrhiz, function(x) x / sum(x))
psrhiz.phylum <- tax_glom(psrhiz.rab, taxrank = "Phylum", NArm = FALSE)

rhiz.phylum <- psmelt(psrhiz.phylum)
rhiz.phylum$Phylum <- as.character(rhiz.phylum$Phylum)
max.rhiz.phylum <- ddply(rhiz.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.phylum <- max.rhiz.phylum[max.rhiz.phylum$max <= .01,]$Phylum
rhiz.phylum[rhiz.phylum$Phylum %in% maxremainder.rhiz.phylum,]$Phylum <- 'Other'

ggplot(rhiz.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Rhizosphere Classes present with at least 1% abundance
```{r rhizclass, echo=FALSE}
psrhiz.class <- tax_glom(psrhiz.rab, taxrank = "Class", NArm = FALSE)

rhiz.class <- psmelt(psrhiz.class)
rhiz.class$Class <- as.character(rhiz.class$Class)
max.rhiz.class <- ddply(rhiz.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.class <- max.rhiz.class[max.rhiz.class$max <= .01,]$Class
rhiz.class[rhiz.class$Class %in% maxremainder.rhiz.class,]$Class <- 'Other'

ggplot(rhiz.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Rhizosphere Families present with at least 2% abundance
```{r rhizfamily, echo=FALSE}
psrhiz.family <- tax_glom(psrhiz.rab, taxrank = "Family", NArm = FALSE)

rhiz.family <- psmelt(psrhiz.family)
rhiz.family$Family <- as.character(rhiz.family$Family)
max.rhiz.family <- ddply(rhiz.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.family <- max.rhiz.family[max.rhiz.family$max <= .02,]$Family
rhiz.family[rhiz.family$Family %in% maxremainder.rhiz.family,]$Family <- 'Other'

ggplot(rhiz.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Rhizosphere genera present with at least 2% abundance
```{r rhizgenera, echo=FALSE}
psrhiz.genus <- tax_glom(psrhiz.rab, taxrank = "Genus", NArm = FALSE)

rhiz.genus <- psmelt(psrhiz.genus)
rhiz.genus$Genus <- as.character(rhiz.genus$Genus)
max.rhiz.genus <- ddply(rhiz.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.genus <- max.rhiz.genus[max.rhiz.genus$max <= .02,]$Genus
rhiz.genus[rhiz.genus$Genus %in% maxremainder.rhiz.genus,]$Genus <- 'Other'

ggplot(rhiz.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()
) + guides(fill=guide_legend(title='Genus')   ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk mid soil Phyla present with at least 1% abundance

```{r echo = FALSE}
psbulkmid.rab <- transform_sample_counts(psbulkmid, function(x) x / sum(x))
psbulkmid.phylum <- tax_glom(psbulkmid.rab, taxrank = "Phylum", NArm = FALSE)

bulkmid.phylum <- psmelt(psbulkmid.phylum)
bulkmid.phylum$Phylum <- as.character(bulkmid.phylum$Phylum)
max.bulkmid.phylum <- ddply(bulkmid.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.phylum <- max.bulkmid.phylum[max.bulkmid.phylum$max <= .01,]$Phylum
bulkmid.phylum[bulkmid.phylum$Phylum %in% maxremainder.bulkmid.phylum,]$Phylum <- 'Other'

ggplot(bulkmid.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), panel.background = element_blank()
) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk mid soil Classes present with at least 2% abundance

```{r echo = FALSE}
psbulkmid.class <- tax_glom(psbulkmid.rab, taxrank = "Class", NArm = FALSE)

bulkmid.class <- psmelt(psbulkmid.class)
bulkmid.class$Class <- as.character(bulkmid.class$Class)
max.bulkmid.class <- ddply(bulkmid.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.class <- max.bulkmid.class[max.bulkmid.class$max <= .02,]$Class
bulkmid.class[bulkmid.class$Class %in% maxremainder.bulkmid.class,]$Class <- 'Other'

ggplot(bulkmid.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk mid soil Families present with at least 3% abundance

```{r echo = FALSE}
psbulkmid.family <- tax_glom(psbulkmid.rab, taxrank = "Family", NArm = FALSE)

bulkmid.family <- psmelt(psbulkmid.family)
bulkmid.family$Family <- as.character(bulkmid.family$Family)
max.bulkmid.family <- ddply(bulkmid.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.family <- max.bulkmid.family[max.bulkmid.family$max <= .03,]$Family
bulkmid.family[bulkmid.family$Family %in% maxremainder.bulkmid.family,]$Family <- 'Other'

ggplot(bulkmid.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk mid soil Genera present with at least 2% abundance

```{r echo = FALSE}
psbulkmid.genus <- tax_glom(psbulkmid.rab, taxrank = "Genus", NArm = FALSE)

bulkmid.genus <- psmelt(psbulkmid.genus)
bulkmid.genus$Genus <- as.character(bulkmid.genus$Genus)
max.bulkmid.genus <- ddply(bulkmid.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.genus <- max.bulkmid.genus[max.bulkmid.genus$max <= .02,]$Genus
bulkmid.genus[bulkmid.genus$Genus %in% maxremainder.bulkmid.genus,]$Genus <- 'Other'

ggplot(bulkmid.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


### Bulk fall soil Phyla present with at least 1% abundance

```{r echo = FALSE}
psbulkfall.rab <- transform_sample_counts(psbulkfall, function(x) x / sum(x))
psbulkfall.phylum <- tax_glom(psbulkfall.rab, taxrank = "Phylum", NArm = FALSE)

bulkfall.phylum <- psmelt(psbulkfall.phylum)
bulkfall.phylum$Phylum <- as.character(bulkfall.phylum$Phylum)
max.bulkfall.phylum <- ddply(bulkfall.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.phylum <- max.bulkfall.phylum[max.bulkfall.phylum$max <= .01,]$Phylum
bulkfall.phylum[bulkfall.phylum$Phylum %in% maxremainder.bulkfall.phylum,]$Phylum <- 'Other'

ggplot(bulkfall.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk fall soil Classes present with at least 2% abundance

```{r echo = FALSE}
psbulkfall.class <- tax_glom(psbulkfall.rab, taxrank = "Class", NArm = FALSE)

bulkfall.class <- psmelt(psbulkfall.class)
bulkfall.class$Class <- as.character(bulkfall.class$Class)
max.bulkfall.class <- ddply(bulkfall.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.class <- max.bulkfall.class[max.bulkfall.class$max <= .02,]$Class
bulkfall.class[bulkfall.class$Class %in% maxremainder.bulkfall.class,]$Class <- 'Other'

ggplot(bulkfall.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk fall soil Families present with at least 3% abundance

```{r echo = FALSE}
psbulkfall.family <- tax_glom(psbulkfall.rab, taxrank = "Family", NArm = FALSE)

bulkfall.family <- psmelt(psbulkfall.family)
bulkfall.family$Family <- as.character(bulkfall.family$Family)
max.bulkfall.family <- ddply(bulkfall.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.family <- max.bulkfall.family[max.bulkfall.family$max <= .03,]$Family
bulkfall.family[bulkfall.family$Family %in% maxremainder.bulkfall.family,]$Family <- 'Other'

ggplot(bulkfall.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"          ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk fall soil Genera present with at least 2% abundance

```{r echo = FALSE}
psbulkfall.genus <- tax_glom(psbulkfall.rab, taxrank = "Genus", NArm = FALSE)

bulkfall.genus <- psmelt(psbulkfall.genus)
bulkfall.genus$Genus <- as.character(bulkfall.genus$Genus)
max.bulkfall.genus <- ddply(bulkfall.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.genus <- max.bulkfall.genus[max.bulkfall.genus$max <= .02,]$Genus
bulkfall.genus[bulkfall.genus$Genus %in% maxremainder.bulkfall.genus,]$Genus <- 'Other'

ggplot(bulkfall.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"                    ) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


### Cyst Mid season Phyla present with at least 1% abundance in cyst sample
```{r cystmidphyla, echo=FALSE}
pscystmid.rab <- transform_sample_counts(pscystmid, function(x) x/sum(x))
pscystmid.phylum <- tax_glom(pscystmid.rab, taxrank = "Phylum", NArm = FALSE)

cystmid.phylum <- psmelt(pscystmid.phylum)
# convert Phylum to a character vector from a factor because R
cystmid.phylum$Phylum <- as.character(cystmid.phylum$Phylum)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystmid.phylum <- ddply(cystmid.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.phylum <- max.cystmid.phylum[max.cystmid.phylum$max <= .01,]$Phylum
cystmid.phylum[cystmid.phylum$Phylum %in% maxremainder.cystmid.phylum,]$Phylum <- 'Other'
ggplot(cystmid.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Cyst Mid season Classes present with at least 1% abundance in cyst sample
```{r cystmidclass, echo=FALSE}
pscystmid.class <- tax_glom(pscystmid.rab, taxrank = "Class", NArm = FALSE)

cystmid.class <- psmelt(pscystmid.class)
# convert Phylum to a character vector from a factor because R
cystmid.class$Class <- as.character(cystmid.class$Class)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystmid.class <- ddply(cystmid.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.class <- max.cystmid.class[max.cystmid.class$max <= .01,]$Class
cystmid.class[cystmid.class$Class %in% maxremainder.cystmid.class,]$Class <- 'Other'
ggplot(cystmid.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), 
panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Cyst Mid season Families present with at least 5% abundance in cyst sample
```{r cystmidfamily, echo=FALSE}
pscystmid.family <- tax_glom(pscystmid.rab, taxrank = "Family", NArm = FALSE)

cystmid.family <- psmelt(pscystmid.family)
# convert taxa to a character vector from a factor because R
cystmid.family$Family <- as.character(cystmid.family$Family)
#Any taxa that is present in a sample greater than 1%
max.cystmid.family <- ddply(cystmid.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.family <- max.cystmid.family[max.cystmid.family$max <= .05,]$Family
cystmid.family[cystmid.family$Family %in% maxremainder.cystmid.family,]$Family <- 'Other'
ggplot(cystmid.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Mid season Genera present with at least 5% abundance in cyst sample
```{r cystmidgenera, echo=FALSE}
pscystmid.genus <- tax_glom(pscystmid.rab, taxrank = "Genus", NArm = FALSE)

cystmid.genus <- psmelt(pscystmid.genus)
# convert taxa to a character vector from a factor because R
cystmid.genus$Genus <- as.character(cystmid.genus$Genus)
#Any taxa that is present in a sample greater than 1%
max.cystmid.genus <- ddply(cystmid.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.genus <- max.cystmid.genus[max.cystmid.genus$max <= .05,]$Genus
cystmid.genus[cystmid.genus$Genus %in% maxremainder.cystmid.genus,]$Genus <- 'Other'

ggplot(cystmid.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')
) +  scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Fall season Phyla present with at least 1% abundance in cyst sample
```{r cystfallphyla, echo=FALSE}
pscystfall.rab <- transform_sample_counts(pscystfall, function(x) x/sum(x))

pscystfall.phylum <- tax_glom(pscystfall.rab, taxrank = "Phylum", NArm = FALSE)

cystfall.phylum <- psmelt(pscystfall.phylum)
# convert Phylum to a character vector from a factor because R
cystfall.phylum$Phylum <- as.character(cystfall.phylum$Phylum)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystfall.phylum <- ddply(cystfall.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.phylum <- max.cystfall.phylum[max.cystfall.phylum$max <= .01,]$Phylum
cystfall.phylum[cystfall.phylum$Phylum %in% maxremainder.cystfall.phylum,]$Phylum <- 'Other'
ggplot(cystfall.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Cyst Fall season Classes present with at least 1% abundance in cyst sample
```{r cystfallclass, echo=FALSE}
pscystfall.class <- tax_glom(pscystfall.rab, taxrank = "Class", NArm = FALSE)

cystfall.class <- psmelt(pscystfall.class)
# convert Phylum to a character vector from a factor because R
cystfall.class$Class <- as.character(cystfall.class$Class)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystfall.class <- ddply(cystfall.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.class <- max.cystfall.class[max.cystfall.class$max <= .01,]$Class
cystfall.class[cystfall.class$Class %in% maxremainder.cystfall.class,]$Class <- 'Other'
ggplot(cystfall.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), 
panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Cyst Fall season Families present with at least 5% abundance in cyst sample
```{r cystfallfamily, echo=FALSE}
pscystfall.family <- tax_glom(pscystfall.rab, taxrank = "Family", NArm = FALSE)

cystfall.family <- psmelt(pscystfall.family)
# convert taxa to a character vector from a factor because R
cystfall.family$Family <- as.character(cystfall.family$Family)
#Any taxa that is present in a sample greater than 1%
max.cystfall.family <- ddply(cystfall.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.family <- max.cystfall.family[max.cystfall.family$max <= .05,]$Family
cystfall.family[cystfall.family$Family %in% maxremainder.cystfall.family,]$Family <- 'Other'
ggplot(cystfall.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Fall season Genera present with at least 5% abundance in cyst sample
```{r cystfallgenus, echo=FALSE}
pscystfall.genus <- tax_glom(pscystfall.rab, taxrank = "Genus", NArm = FALSE)

cystfall.genus <- psmelt(pscystfall.genus)
# convert taxa to a character vector from a factor because R
cystfall.genus$Genus <- as.character(cystfall.genus$Genus)
#Any taxa that is present in a sample greater than 1%
max.cystfall.genus <- ddply(cystfall.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.genus <- max.cystfall.genus[max.cystfall.genus$max <= .05,]$Genus
cystfall.genus[cystfall.genus$Genus %in% maxremainder.cystfall.genus,]$Genus <- 'Other'

ggplot(cystfall.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')
) +  scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Root Phyla present in the root samples
```{r rootphyla, echo=FALSE}
psroot.rab <- transform_sample_counts(psroot, function(x) x/sum(x))

psroot.phylum <- tax_glom(psroot.rab, taxrank = "Phylum", NArm = FALSE)
root.phylum <- psmelt(psroot.phylum)
# convert Phylum to a character vector from a factor because R
root.phylum$Phylum <- as.character(root.phylum$Phylum)
ggplot(root.phylum, aes(x = Sample, y = Abundance, fill = Phylum))+ geom_bar(stat = "identity"
)+ facet_grid(~Plot, scales = "free", space = "free") + theme(axis.text.x = element_blank(), 
strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Root Classes present with at least a 1% abundance in root samples
```{r rootclass, echo=FALSE}
psroot.class <- tax_glom(psroot.rab, taxrank = "Class", NArm = FALSE)

root.class <- psmelt(psroot.class)
# convert Phylum to a character vector from a factor because R
root.class$Class <- as.character(root.class$Class)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.root.class <- ddply(root.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.root.class <- max.root.class[max.root.class$max <= .01,]$Class
root.class[root.class$Class %in% maxremainder.root.class,]$Class <- 'Other'

ggplot(root.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Root Families present with at least 3% abundance in root samples
```{r rootfamily, echo=FALSE}
psroot.family <- tax_glom(psroot.rab, taxrank = "Family", NArm = FALSE)

root.family <- psmelt(psroot.family)
# convert Phylum to a character vector from a factor because R
root.family$Family <- as.character(root.family$Family)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.root.family <- ddply(root.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.root.family <- max.root.family[max.root.family$max <= .03,]$Family
root.family[root.family$Family %in% maxremainder.root.family,]$Family <- 'Other'

ggplot(root.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Root Genera present with at least 3% abundance in root samples
```{r rootgenus, echo=FALSE}
psroot.genus <- tax_glom(psroot.rab, taxrank = "Genus", NArm=FALSE)

# create dataframe from phyloseq object
root.genus <- psmelt(psroot.genus)
# convert Phylum to a character vector from a factor because R
root.genus$Genus <- as.character(root.genus$Genus)
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.root.genus <- ddply(root.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.root.genus <- max.root.genus[max.root.genus$max <= .03,]$Genus
root.genus[root.genus$Genus %in% maxremainder.root.genus,]$Genus <- 'Other'

ggplot(root.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

## Investigating the "Core" Microbiome 

### Entire dataset. Taxa present at .01% abundance or higher in 50% of the samples
```{r  echo=FALSE}
psall <- ps4
psall.rab <- transform_sample_counts(psall, function(x) x/sum(x))
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .50*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Entire data set core microbiome. Genera present with at least .01% abundance in 60% of the samples
```{r  echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .6*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 10)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Entire data set core microbiome. Genera present with at least .01% abundance in  70% of samples
```{r  echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .70*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(
),strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

```

### Entire data set core microbiome. Genera present with at least .01% abundance in 75% of the samples
```{r echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .75*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(
),strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Entire data set core microbiome. Genera present with at least .01% abundance in 80% of the samples
```{r echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .80*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(
),strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Entire data set core microbiome. Genera present with at least .01% abundance in 90% of the samples
```{r echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .90*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(
),strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Split rhizosphere into plots
```{r}
psrhizSa <- subset_samples(psrhiz, Plot == "Sa")
psrhizS3 <- subset_samples(psrhiz, Plot == "S3")
psrhizSs <- subset_samples(psrhiz, Plot == "Ss")
psrhizSr <- subset_samples(psrhiz, Plot == "Sr")
```

### SA Rhizosphere (1% abundance or higher) 
```{r  echo=FALSE}
psrhizSa.rab <- transform_sample_counts(psrhizSa, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizSa.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizSacore.rab <- prune_taxa(filter.rab, psrhizSa.rab)
rhizSacore.genus.rab <- tax_glom(rhizSacore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizSacore.genus.rab.df <- psmelt(rhizSacore.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSacore.genus.rab.df$Genus <- as.character(rhizSacore.genus.rab.df$Genus)
ggplot(rhizSacore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### S3 Rhizosphere (1% abundance or higher)
```{r  echo=FALSE}
psrhizS3.rab <- transform_sample_counts(psrhizS3, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizS3.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizS3core.rab <- prune_taxa(filter.rab, psrhizS3.rab)
rhizS3core.genus.rab <- tax_glom(rhizS3core.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizS3core.genus.rab.df <- psmelt(rhizS3core.genus.rab)
# convert Genus to a character vector from a factor because R
rhizS3core.genus.rab.df$Genus <- as.character(rhizS3core.genus.rab.df$Genus)
ggplot(rhizS3core.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### SR Rhizosphere (1% abundance or higher)
```{r  echo=FALSE}
psrhizSr.rab <- transform_sample_counts(psrhizSr, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizSr.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizSrcore.rab <- prune_taxa(filter.rab, psrhizSr.rab)
rhizSrcore.genus.rab <- tax_glom(rhizSrcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizSrcore.genus.rab.df <- psmelt(rhizSrcore.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSrcore.genus.rab.df$Genus <- as.character(rhizSrcore.genus.rab.df$Genus)
ggplot(rhizSrcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### SS Rhizosphere (1% abundance or higher)
```{r  echo=FALSE}
psrhizSs.rab <- transform_sample_counts(psrhizSs, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizSs.rab, filterfun_sample(function(x) x / sum(x) > .01))
rhizSscore.rab <- prune_taxa(filter.rab, psrhizSs.rab)
rhizSscore.genus.rab <- tax_glom(rhizSscore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizSscore.genus.rab.df <- psmelt(rhizSscore.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSscore.genus.rab.df$Genus <- as.character(rhizSscore.genus.rab.df$Genus)
ggplot(rhizSscore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Split bulk soil samples into plots & seasons
```{r  echo=FALSE}
psbulkSa <- subset_samples(psbulk, Plot == "Sa")
psbulkS3 <- subset_samples(psbulk, Plot == "S3")
psbulkSs <- subset_samples(psbulk, Plot == "Ss")
psbulkSr <- subset_samples(psbulk, Plot == "Sr")

#split into mid season
psbulkSamid <- subset_samples(psbulkSa, Season == "Mid")
psbulkS3mid <- subset_samples(psbulkS3, Season == "Mid")
psbulkSsmid <- subset_samples(psbulkSs, Season == "Mid")
psbulkSrmid <- subset_samples(psbulkSr, Season == "Mid")

psbulkSafall <- subset_samples(psbulkSa, Season == "Fall")
psbulkS3fall <- subset_samples(psbulkS3, Season == "Fall")
psbulkSsfall <- subset_samples(psbulkSs, Season == "Fall")
psbulkSrfall <- subset_samples(psbulkSr, Season == "Fall")
```

### Bulk soil Mid season SA core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkSamid.rab <- transform_sample_counts(psbulkSamid, function(x) x/sum(x))
filterbSamid.rab <- phyloseq::genefilter_sample(psbulkSamid.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSamidcore.rab <- prune_taxa(filterbSamid.rab, psbulkSamid.rab)
bulkSamidcore.genus.rab <- tax_glom(bulkSamidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSamidcore.genus.rab.df <- psmelt(bulkSamidcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSamidcore.genus.rab.df$Genus <- as.character(bulkSamidcore.genus.rab.df$Genus)
ggplot(bulkSamidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk soil Mid season S3 core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkS3mid.rab <- transform_sample_counts(psbulkS3mid, function(x) x/sum(x))
filterbS3mid.rab <- phyloseq::genefilter_sample(psbulkS3mid.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkS3midcore.rab <- prune_taxa(filterbS3mid.rab, psbulkS3mid.rab)
bulkS3midcore.genus.rab <- tax_glom(bulkS3midcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkS3midcore.genus.rab.df <- psmelt(bulkS3midcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3midcore.genus.rab.df$Genus <- as.character(bulkS3midcore.genus.rab.df$Genus)
ggplot(bulkS3midcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk soil Mid season SR core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkSrmid.rab <- transform_sample_counts(psbulkSrmid, function(x) x/sum(x))
filterbSrmid.rab <- phyloseq::genefilter_sample(psbulkSrmid.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSrmidcore.rab <- prune_taxa(filterbSrmid.rab, psbulkSrmid.rab)
bulkSrmidcore.genus.rab <- tax_glom(bulkSrmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSrmidcore.genus.rab.df <- psmelt(bulkSrmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSrmidcore.genus.rab.df$Genus <- as.character(bulkSrmidcore.genus.rab.df$Genus)
ggplot(bulkSrmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk soil Mid season SS core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkSsmid.rab <- transform_sample_counts(psbulkSsmid, function(x) x/sum(x))
filterbSsmid.rab <- phyloseq::genefilter_sample(psbulkSsmid.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSsmidcore.rab <- prune_taxa(filterbSsmid.rab, psbulkSsmid.rab)
bulkSsmidcore.genus.rab <- tax_glom(bulkSsmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSsmidcore.genus.rab.df <- psmelt(bulkSsmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSsmidcore.genus.rab.df$Genus <- as.character(bulkSsmidcore.genus.rab.df$Genus)
ggplot(bulkSsmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk soil SA Fall season core genera (1% abundance or higher) 
```{r  echo=FALSE}
psbulkSafall.rab <- transform_sample_counts(psbulkSafall, function(x) x/sum(x))
filterbSafall.rab <- phyloseq::genefilter_sample(psbulkSafall.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSafallcore.rab <- prune_taxa(filterbSafall.rab, psbulkSafall.rab)
bulkSafallcore.genus.rab <- tax_glom(bulkSafallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSafallcore.genus.rab.df <- psmelt(bulkSafallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSafallcore.genus.rab.df$Genus <- as.character(bulkSafallcore.genus.rab.df$Genus)
ggplot(bulkSafallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk soil S3 fall season core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkS3fall.rab <- transform_sample_counts(psbulkS3fall, function(x) x/sum(x))
filterbS3fall.rab <- phyloseq::genefilter_sample(psbulkS3fall.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkS3fallcore.rab <- prune_taxa(filterbS3fall.rab, psbulkS3fall.rab)
bulkS3fallcore.genus.rab <- tax_glom(bulkS3fallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkS3fallcore.genus.rab.df <- psmelt(bulkS3fallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3fallcore.genus.rab.df$Genus <- as.character(bulkS3fallcore.genus.rab.df$Genus)
ggplot(bulkS3fallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


### Bulk soil SR fall season core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkSrfall.rab <- transform_sample_counts(psbulkSrfall, function(x) x/sum(x))
filterbSrfall.rab <- phyloseq::genefilter_sample(psbulkSrfall.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSrfallcore.rab <- prune_taxa(filterbSrfall.rab, psbulkSrfall.rab)
bulkSrfallcore.genus.rab <- tax_glom(bulkSrfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSrfallcore.genus.rab.df <- psmelt(bulkSrfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSrfallcore.genus.rab.df$Genus <- as.character(bulkSrfallcore.genus.rab.df$Genus)
ggplot(bulkSrfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk soil SS fall season core genera (1% abundance or higher)
```{r  echo=FALSE}
psbulkSsfall.rab <- transform_sample_counts(psbulkSsfall, function(x) x/sum(x))
filterbSsfall.rab <- phyloseq::genefilter_sample(psbulkSsfall.rab, filterfun_sample(function(x) x / sum(x) > .01))
bulkSsfallcore.rab <- prune_taxa(filterbSsfall.rab, psbulkSsfall.rab)
bulkSsfallcore.genus.rab <- tax_glom(bulkSsfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSsfallcore.genus.rab.df <- psmelt(bulkSsfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSsfallcore.genus.rab.df$Genus <- as.character(bulkSsfallcore.genus.rab.df$Genus)
ggplot(bulkSsfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Split cyst samples in plots and seasons

```{r  echo=FALSE}
pscystSa <- subset_samples(pscyst, Plot == "Sa")
pscystS3 <- subset_samples(pscyst, Plot == "S3")
pscystSs <- subset_samples(pscyst, Plot == "Ss")
pscystSr <- subset_samples(pscyst, Plot == "Sr")

#split into midseason
pscystSamid <- subset_samples(pscystSa, Season == "Mid")
pscystS3mid <- subset_samples(pscystS3, Season == "Mid")
pscystSsmid <- subset_samples(pscystSs, Season == "Mid")
pscystSrmid <- subset_samples(pscystSr, Season == "Mid")

pscystSafall <- subset_samples(pscystSa, Season == "Fall")
pscystS3fall <- subset_samples(pscystS3, Season == "Fall")
pscystSsfall <- subset_samples(pscystSs, Season == "Fall")
pscystSrfall <- subset_samples(pscystSr, Season == "Fall")
```

### Cysts SA mid season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystSamid.rab <- transform_sample_counts(pscystSamid, function(x) x/sum(x))
filtercSamid.rab <- phyloseq::genefilter_sample(pscystSamid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSamidcore.rab <- prune_taxa(filtercSamid.rab, pscystSamid.rab)
cystSamidcore.genus.rab <- tax_glom(cystSamidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSamidcore.genus.rab.df <- psmelt(cystSamidcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSamidcore.genus.rab.df$Genus <- as.character(cystSamidcore.genus.rab.df$Genus)
ggplot(cystSamidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts S3 mid season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystS3mid.rab <- transform_sample_counts(pscystS3mid, function(x) x/sum(x))
filtercS3mid.rab <- phyloseq::genefilter_sample(pscystS3mid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3midcore.rab <- prune_taxa(filtercS3mid.rab, pscystS3mid.rab)
cystS3midcore.genus.rab <- tax_glom(cystS3midcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystS3midcore.genus.rab.df <- psmelt(cystS3midcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3midcore.genus.rab.df$Genus <- as.character(cystS3midcore.genus.rab.df$Genus)
ggplot(cystS3midcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts SR mid season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystSrmid.rab <- transform_sample_counts(pscystSrmid, function(x) x/sum(x))
filtercSrmid.rab <- phyloseq::genefilter_sample(pscystSrmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSrmidcore.rab <- prune_taxa(filtercSrmid.rab, pscystSrmid.rab)
cystSrmidcore.genus.rab <- tax_glom(cystSrmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSrmidcore.genus.rab.df <- psmelt(cystSrmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSrmidcore.genus.rab.df$Genus <- as.character(cystSrmidcore.genus.rab.df$Genus)
ggplot(cystSrmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts SS mid season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystSsmid.rab <- transform_sample_counts(pscystSsmid, function(x) x/sum(x))
filtercSsmid.rab <- phyloseq::genefilter_sample(pscystSsmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSsmidcore.rab <- prune_taxa(filtercSsmid.rab, pscystSsmid.rab)
cystSsmidcore.genus.rab <- tax_glom(cystSsmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSsmidcore.genus.rab.df <- psmelt(cystSsmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSsmidcore.genus.rab.df$Genus <- as.character(cystSsmidcore.genus.rab.df$Genus)
ggplot(cystSsmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts SA fall season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystSafall.rab <- transform_sample_counts(pscystSafall, function(x) x/sum(x))
filtercSafall.rab <- phyloseq::genefilter_sample(pscystSafall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSafallcore.rab <- prune_taxa(filtercSafall.rab, pscystSafall.rab)
cystSafallcore.genus.rab <- tax_glom(cystSafallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSafallcore.genus.rab.df <- psmelt(cystSafallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSafallcore.genus.rab.df$Genus <- as.character(cystSafallcore.genus.rab.df$Genus)
ggplot(cystSafallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

```

### Cysts S3 fall season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystS3fall.rab <- transform_sample_counts(pscystS3fall, function(x) x/sum(x))
filtercS3fall.rab <- phyloseq::genefilter_sample(pscystS3fall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3fallcore.rab <- prune_taxa(filtercS3fall.rab, pscystS3fall.rab)
cystS3fallcore.genus.rab <- tax_glom(cystS3fallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystS3fallcore.genus.rab.df <- psmelt(cystS3fallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3fallcore.genus.rab.df$Genus <- as.character(cystS3fallcore.genus.rab.df$Genus)
ggplot(cystS3fallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts SR fall season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystSrfall.rab <- transform_sample_counts(pscystSrfall, function(x) x/sum(x))
filtercSrfall.rab <- phyloseq::genefilter_sample(pscystSrfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSrfallcore.rab <- prune_taxa(filtercSrfall.rab, pscystSrfall.rab)
cystSrfallcore.genus.rab <- tax_glom(cystSrfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSrfallcore.genus.rab.df <- psmelt(cystSrfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSrfallcore.genus.rab.df$Genus <- as.character(cystSrfallcore.genus.rab.df$Genus)
ggplot(cystSrfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts SS fall season core genera (2% abundance or higher)
```{r  echo=FALSE}
pscystSsfall.rab <- transform_sample_counts(pscystSsfall, function(x) x/sum(x))
filtercSsfall.rab <- phyloseq::genefilter_sample(pscystSsfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSsfallcore.rab <- prune_taxa(filtercSsfall.rab, pscystSsfall.rab)
cystSsfallcore.genus.rab <- tax_glom(cystSsfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSsfallcore.genus.rab.df <- psmelt(cystSsfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSsfallcore.genus.rab.df$Genus <- as.character(cystSsfallcore.genus.rab.df$Genus)
ggplot(cystSsfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Split roots into plots
```{r  echo=FALSE}
psrootSa <- subset_samples(psroot, Plot == "Sa")
psrootS3 <- subset_samples(psroot, Plot == "S3")
psrootSs <- subset_samples(psroot, Plot == "Ss")
psrootSr <- subset_samples(psroot, Plot == "Sr")
```

### Roots SA core genera (2% abundance or higher)
```{r  echo=FALSE}
psrootSa.rab <- transform_sample_counts(psrootSa, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootSa.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootSacore.rab <- prune_taxa(filter.rab, psrootSa.rab)
rootSacore.genus.rab <- tax_glom(rootSacore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootSacore.genus.rab.df <- psmelt(rootSacore.genus.rab)
# convert Genus to a character vector from a factor because R
rootSacore.genus.rab.df$Genus <- as.character(rootSacore.genus.rab.df$Genus)
ggplot(rootSacore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Roots S3 core genera (2% abundance or higher)
```{r  echo=FALSE}
psrootS3.rab <- transform_sample_counts(psrootS3, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootS3.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootS3core.rab <- prune_taxa(filter.rab, psrootS3.rab)
rootS3core.genus.rab <- tax_glom(rootS3core.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootS3core.genus.rab.df <- psmelt(rootS3core.genus.rab)
# convert Genus to a character vector from a factor because R
rootS3core.genus.rab.df$Genus <- as.character(rootS3core.genus.rab.df$Genus)
ggplot(rootS3core.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Roots SR core genera (2% abundance or higher)
```{r  echo=FALSE}
psrootSr.rab <- transform_sample_counts(psrootSr, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootSr.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootSrcore.rab <- prune_taxa(filter.rab, psrootSr.rab)
rootSrcore.genus.rab <- tax_glom(rootSrcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootSrcore.genus.rab.df <- psmelt(rootSrcore.genus.rab)
# convert Genus to a character vector from a factor because R
rootSrcore.genus.rab.df$Genus <- as.character(rootSrcore.genus.rab.df$Genus)
ggplot(rootSrcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Roots SS core genera (2% abundance or higher) 
```{r  echo=FALSE}
psrootSs.rab <- transform_sample_counts(psrootSs, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootSs.rab, filterfun_sample(function(x) x / sum(x) > .02))
rootSscore.rab <- prune_taxa(filter.rab, psrootSs.rab)
rootSscore.genus.rab <- tax_glom(rootSscore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootSscore.genus.rab.df <- psmelt(rootSscore.genus.rab)
# convert Genus to a character vector from a factor because R
rootSscore.genus.rab.df$Genus <- as.character(rootSscore.genus.rab.df$Genus)
ggplot(rootSscore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


# SCN ITS
---
title: "LTR_ITS_S3ARSMarkdown"
author: "Emily Green"
date: "2023-11-03"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Upload data and assemble reads

### Run DADA2

```{r dada2, echo = FALSE, results = "hide", include=FALSE, message=FALSE, cache=TRUE}
library(dada2)
path <- "/workdir/eag252/LTR_ITS_S3ARS/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
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
#It looks like a lot of short reads were filtered out. The minLen was 50. Why would I have so many short reads?

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
#write.csv(track, "allITScounts.csv")


# UNITE fungal database from https://doi.plutof.ut.ee/doi/10.15156/BIO/2483911
# sh_general_release_27.10.2022.tgz
taxa <- assignTaxonomy(seqtab.nochim, "/workdir/eag252/LTR_ITS_S3ARS/sh_general_release_dynamic_27.10.2022.fasta", multithread=TRUE)

#taxa.print <- taxa 
#rownames(taxa.print) <- NULL
#head(taxa.print)

```

### Load packages
```{r loadpackages, message=FALSE, echo=FALSE}
library(phyloseq)
library(ggplot2)
library(plyr)
library(vegan)
library(viridisLite)
library(viridis)
library(ape)
library(cowplot)
theme_set(theme_bw())
```

## Import metadata
```{r importmetadata, echo=FALSE}
meta <- read.csv("/workdir/eag252/LTR_ITS_S3ARS/LTRITS_S3ARSMetadata.csv",header = TRUE, row.names = 1)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(meta), tax_table(taxa))
#sample_data(ps)

```

## Quality Control

### Filter taxa

```{r filttaxa, echo=FALSE, results = "hide"}
sort(get_taxa_unique(ps, taxonomic.rank = "Kingdom"))
sort(get_taxa_unique(ps, taxonomic.rank = "Phylum"))
#View(tax_table(ps))
ps1 <- subset_taxa(ps, !Phylum == "NA")

#View(tax_table(ps1))
```

### Run Decontam
```{r decontam, echo=FALSE}
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

#sample_sums(ps1)
#View(sample_sums(ps1))
```

### Quick check of alpha diversity
```{r quickalpha, echo=FALSE}
plot_richness(ps1, x="SampleType", measures=c("Shannon", "Simpson"), color="Sample_or_Control")
```

### Remove samples with less than 100 reads
```{r filtsamp, echo=FALSE}
ps2 <- prune_samples(sample_sums(ps1) > 100, ps1)
```

### Beta Diversity Check if controls and samples cluster separately
```{r betacheck, echo=FALSE}
ps2t <- ps2
#Unifrac
ps2t.tree <- rtree(ntaxa(ps2t), rooted = TRUE, tip.label = taxa_names(ps2t))
#plot(ps2t.tree)
ps2t <- merge_phyloseq(ps2t, ps2t.tree)
```

#### Unweighted

```{r unweight2, echo=FALSE}
ps2tunifrac_dist <- phyloseq::distance(ps2t, method="unifrac", weighted=F)
ordination = ordinate(ps2t, method="PCoA", distance=ps2tunifrac_dist)
plot_ordination(ps2t, ordination, color = "Sample_or_Control", shape = "SampleType") + theme(aspect.ratio=1)
```

```{r}
ps2t.tree.data <- data.frame(sample_data(ps2t))
adonis2(ps2tunifrac_dist ~SampleType*Sample_or_Control, ps2t.tree.data)
```

#### Weighted
```{r weight2, echo=FALSE}
ps2tunifrac_dist <- phyloseq::distance(ps2t, method="unifrac", weighted=T)
ordination = ordinate(ps2t, method="PCoA", distance=ps2tunifrac_dist)
plot_ordination(ps2t, ordination, color = "Sample_or_Control", shape = "SampleType") + theme(aspect.ratio=1)
```

```{r}
ps2t.tree.data <- data.frame(sample_data(ps2t))
adonis2(ps2tunifrac_dist ~SampleType*Sample_or_Control, ps2t.tree.data)
```

### Remove all negative controls and mock communities
```{r filtsample, echo=FALSE}
ps3 <- subset_samples(ps2, Sample_or_Control == "Sample")
```

### Quick view of beta-diversity of soil, cysts and roots.
```{r betaofall, echo=FALSE}
ps3t <- ps3
ps3t.tree <- rtree(ntaxa(ps3t), rooted = TRUE, tip.label = taxa_names(ps3t))
#plot(ps3t.tree)
ps3t <- merge_phyloseq(ps3t, ps3t.tree)
```

#### Unweighted
```{r unweight3, echo=FALSE}
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=F)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "SampleType", shape = "Season") + theme(aspect.ratio=1)
```

```{r}
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~SampleType*Season, ps3t.tree.data)
```

#### Weighted
```{r weight3, echo=FALSE}
ps3tunifrac_dist <- phyloseq::distance(ps3t, method="unifrac", weighted=T)
ordination = ordinate(ps3t, method="PCoA", distance=ps3tunifrac_dist)
plot_ordination(ps3t, ordination, color = "SampleType", shape = "Season") + theme(aspect.ratio=1)
```

```{r}
ps3t.tree.data <- data.frame(sample_data(ps3t))
adonis2(ps3tunifrac_dist ~SampleType*Season, ps3t.tree.data)
```

## Split samples based on 'Sample Type' into separate phyloseq objects
```{r splitsamples}
pssoil <- subset_samples(ps3, SampleType == "Soil")
pscyst <- subset_samples(ps3, SampleType == "Cyst")
psroot <- subset_samples(ps3, SampleType == "Roots")
```

### Confirm that each phyloseq object contains the correct samples
```{r checkpsobject }
#View(sample_data(pssoil))
#View(sample_data(pscyst))
#View(sample_data(psroot))
```

### Remove taxa with 0 reads in each phyloseq object
```{r removereads}
pssoil <- prune_taxa(taxa_sums(pssoil) > 0, pssoil)
pscyst <- prune_taxa(taxa_sums(pscyst) > 0, pscyst)
psroot <- prune_taxa(taxa_sums(psroot) > 0, psroot)
```

## Alpha Diversity & Anova statistics

### Soil Alpha Diversity

```{r soilalpha, echo=FALSE}
soilalpha <- estimate_richness(pssoil, measures = c("Shannon", "Simpson"))
soilalpha$Plot <- sample_data(pssoil)$Plot
soilalpha$SoilType <- sample_data(pssoil)$SoilType
soilalpha$Season <- sample_data(pssoil)$Season
```

#### Soil Shannon Alpha Diversity
```{r echo=FALSE}
ggplot(soilalpha, aes(x = Plot, y = Shannon, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
```

```{r soilshan}
soilalpha.soiltype <- aov(Shannon ~SoilType, soilalpha)
anova(soilalpha.soiltype)

soilalpha.season <- aov(Shannon ~Season, soilalpha)
anova(soilalpha.season)

soilalpha.plot <- aov(Shannon ~Plot, soilalpha)
anova(soilalpha.plot)
```

#### Soil Simpson Alpha Diversity
```{r echo=FALSE}
ggplot(soilalpha, aes(x = Plot, y = Simpson, color = Season)) + geom_boxplot()+ facet_grid(~SoilType, scales = "free", space = "free" )
```

```{r soilsimp}
soilalpha.soiltype <- aov(Simpson ~SoilType, soilalpha)
anova(soilalpha.soiltype)

soilalpha.season <- aov(Simpson ~Season, soilalpha)
anova(soilalpha.season)

soilalpha.plot <- aov(Simpson ~Plot, soilalpha)
anova(soilalpha.plot)
```

### Cyst Alpha Diversity
```{r cystalpha, echo=FALSE}
cystalpha <- estimate_richness(pscyst,measures = c("Shannon", "Simpson"))
cystalpha$Plot <- sample_data(pscyst)$Plot
cystalpha$Season <- sample_data(pscyst)$Season
```

#### Cyst Shannon Alpha Diversity
```{r echo=FALSE}
ggplot(cystalpha, aes(x= Plot, y = Shannon, color = Season)) + geom_boxplot()
```

```{r cystshan}
cystalpha.plot <- aov(Shannon ~Plot, cystalpha)
anova(cystalpha.plot)

cystalpha.season <- aov(Shannon ~Season, cystalpha)
anova(cystalpha.season)
```

#### Cyst Simpson Alpha Diversity
```{r echo=FALSE}
ggplot(cystalpha, aes(x= Plot, y = Simpson, color = Season)) + geom_boxplot()
```

```{r cystsimp}
cystalpha.plot <- aov(Simpson ~Plot, cystalpha)
anova(cystalpha.plot)
cystalpha.season <- aov(Simpson ~Season, cystalpha)
anova(cystalpha.season)
```


### Root Alpha Diversity
```{r rootalpha, echo=FALSE}
rootalpha <- estimate_richness(psroot, measures = c("Shannon", "Simpson"))
rootalpha$Plot <- sample_data(psroot)$Plot
```

#### Root Shannon Alpha Diversity
```{r echo=FALSE}
ggplot(rootalpha, aes(x= Plot, y = Shannon)) + geom_boxplot()
```

```{r rootshan}
rootalpha.plot <- aov(Shannon ~Plot, rootalpha)
anova(rootalpha.plot)
```

#### Root Simpson Alpha Diversity
```{r echo=FALSE}
ggplot(rootalpha, aes(x= Plot, y = Simpson)) + geom_boxplot()
```

```{r rootsimp}
rootalpha.plot <- aov(Simpson ~Plot, rootalpha)
anova(rootalpha.plot)
```

## Beta Diversity & PERMANOVA


#### Create phyloseq files for unifrac & split samples into seasons
```{r subsetunifracps}
pssoil.t <- pssoil
pscyst.t <- pscyst
psroot.t <- psroot

psrhiz <- subset_samples(pssoil, SoilType == "Rhizosphere")
psrhiz.t <- psrhiz
psbulk <- subset_samples(pssoil, SoilType == "Bulk soil")
psbulk.t <- psbulk

psbulkmid <- subset_samples(psbulk, Season == "Mid")
psbulkmid.t <- psbulkmid
psbulkfall <- subset_samples(psbulk, Season == "Fall")
psbulkfall.t <- psbulkfall

pscystmid <- subset_samples(pscyst, Season == "Mid")
pscystmid.t <- pscystmid
pscystfall <- subset_samples(pscyst, Season == "Fall")
pscystfall.t <- pscystfall
```

### Soil Beta Diversity
```{r soilbeta, echo=FALSE}
pssoil.t.tree <- rtree(ntaxa(pssoil.t), rooted = TRUE, tip.label = taxa_names(pssoil.t))
#plot(pssoil.t.tree)
pssoil.t <- merge_phyloseq(pssoil.t, pssoil.t.tree)
```

#### Soil Unweighted
```{r soilunweight, echo=FALSE}
pssoiltunifrac_dist <- phyloseq::distance(pssoil.t, method="unifrac", weighted=F)
ordination = ordinate(pssoil.t, method="PCoA", distance=pssoiltunifrac_dist)
plot_ordination(pssoil.t, ordination, color = "SoilType", shape = "Season", label = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.t, ordination, color = "Plot", shape = "SoilType", label = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pssoil.t.tree.data <- data.frame(sample_data(pssoil.t))
adonis2(pssoiltunifrac_dist ~SoilType*Plot*Season, pssoil.t.tree.data)
```

#### Soil Weighted

```{r soilweight, echo=FALSE}
pssoiltunifrac_dist <- phyloseq::distance(pssoil.t, method="unifrac", weighted=T)
ordination = ordinate(pssoil.t, method="PCoA", distance=pssoiltunifrac_dist)
plot_ordination(pssoil.t, ordination, color = "SoilType", shape = "Season", label = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pssoil.t, ordination, color = "Plot", shape = "SoilType", label = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pssoil.t.tree.data <- data.frame(sample_data(pssoil.t))
adonis2(pssoiltunifrac_dist ~SoilType*Plot*Season, pssoil.t.tree.data)
```

### Rhizosphere Beta Diversity
```{r rhizbeta, echo=FALSE}
psrhiz.t.tree <- rtree(ntaxa(psrhiz.t), rooted = TRUE, tip.label = taxa_names(psrhiz.t))
#plot(psrhiz.t.tree)
psrhiz.t <- merge_phyloseq(psrhiz.t, psrhiz.t.tree)
```

#### Rhizosphere Unweighted
```{r rhizunweight, echo=FALSE}
psrhiztunifrac_dist <- phyloseq::distance(psrhiz.t, method="unifrac", weighted=F)
ordination = ordinate(psrhiz.t, method="PCoA", distance=psrhiztunifrac_dist)
plot_ordination(psrhiz.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psrhiz.t.tree.data <- data.frame(sample_data(psrhiz.t))
adonis2(psrhiztunifrac_dist ~Plot, psrhiz.t.tree.data)
```

#### Rhizosphere Weighted
```{r rhizweight, echo=FALSE}
psrhiztunifrac_dist <- phyloseq::distance(psrhiz.t, method="unifrac", weighted=T)
ordination = ordinate(psrhiz.t, method="PCoA", distance=psrhiztunifrac_dist)
plot_ordination(psrhiz.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psrhiz.t.tree.data <- data.frame(sample_data(psrhiz.t))
adonis2(psrhiztunifrac_dist ~Plot, psrhiz.t.tree.data)
```

### Bulk Soil Beta Diversity
```{r bulkbeta, echo=FALSE}
psbulk.t.tree <- rtree(ntaxa(psbulk.t), rooted = TRUE, tip.label = taxa_names(psbulk.t))
#plot(psbulk.t.tree)
psbulk.t <- merge_phyloseq(psbulk.t, psbulk.t.tree)
```

#### Bulk Soil Unweighted
```{r bulkunweight, echo=FALSE}
psbulktunifrac_dist <- phyloseq::distance(psbulk.t, method="unifrac", weighted=F)
ordination = ordinate(psbulk.t, method="PCoA", distance=psbulktunifrac_dist)
plot_ordination(psbulk.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulk.t.tree.data <- data.frame(sample_data(psbulk.t))
adonis2(psbulktunifrac_dist ~Plot*Season, psbulk.t.tree.data)
```

#### Bulk Soil Weighted
```{r bulkweight, echo=FALSE}
psbulktunifrac_dist <- phyloseq::distance(psbulk.t, method="unifrac", weighted=T)
ordination = ordinate(psbulk.t, method="PCoA", distance=psbulktunifrac_dist)
plot_ordination(psbulk.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulk.t.tree.data <- data.frame(sample_data(psbulk.t))
adonis2(psbulktunifrac_dist ~Plot*Season, psbulk.t.tree.data)
```

### Bulk soil Mid season Beta Diversity
```{r bulkmidbeta, echo=FALSE}
psbulkmid.t.tree <- rtree(ntaxa(psbulkmid.t), rooted = TRUE, tip.label = taxa_names(psbulkmid.t))
#plot(psbulkmid.t.tree)
psbulkmid.t <- merge_phyloseq(psbulkmid.t, psbulkmid.t.tree)
```

#### Unweighted Bulk Mid
```{r bulkmidunweight, echo=FALSE}
psbulkmidtunifrac_dist <- phyloseq::distance(psbulkmid.t, method="unifrac", weighted=F)
ordination = ordinate(psbulkmid.t, method="PCoA", distance=psbulkmidtunifrac_dist)
plot_ordination(psbulkmid.t, ordination, color = "Plot",) + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkmid.t.tree.data <- data.frame(sample_data(psbulkmid.t))
adonis2(psbulkmidtunifrac_dist ~Plot, psbulkmid.t.tree.data)
```

#### Weighted Bulk Mid
```{r bulkmidweight, echo=FALSE}
psbulkmidtunifrac_dist <- phyloseq::distance(psbulkmid.t, method="unifrac", weighted=T)
ordination = ordinate(psbulkmid.t, method="PCoA", distance=psbulkmidtunifrac_dist)
plot_ordination(psbulkmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkmid.t.tree.data <- data.frame(sample_data(psbulkmid.t))
adonis2(psbulkmidtunifrac_dist ~Plot, psbulkmid.t.tree.data)
```

### Bulk soil Fall season Beta Diversity
```{r bulkfallbeta, echo=FALSE}
psbulkfall.t.tree <- rtree(ntaxa(psbulkfall.t), rooted = TRUE, tip.label = taxa_names(psbulkfall.t))
#plot(psbulkfall.t.tree)
psbulkfall.t <- merge_phyloseq(psbulkfall.t, psbulkfall.t.tree)
```

#### Unweighted Bulk Fall
```{r bulkfallunweight, echo=FALSE}
psbulkfalltunifrac_dist <- phyloseq::distance(psbulkfall.t, method="unifrac", weighted=F)
ordination = ordinate(psbulkfall.t, method="PCoA", distance=psbulkfalltunifrac_dist)
plot_ordination(psbulkfall.t, ordination, color = "Plot",) + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkfall.t.tree.data <- data.frame(sample_data(psbulkfall.t))
adonis2(psbulkfalltunifrac_dist ~Plot, psbulkfall.t.tree.data)
```

#### Weighted Bulk Fall
```{r bulkfallweight, echo=FALSE}
psbulkfalltunifrac_dist <- phyloseq::distance(psbulkfall.t, method="unifrac", weighted=T)
ordination = ordinate(psbulkfall.t, method="PCoA", distance=psbulkfalltunifrac_dist)
plot_ordination(psbulkfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psbulkfall.t.tree.data <- data.frame(sample_data(psbulkfall.t))
adonis2(psbulkfalltunifrac_dist ~Plot, psbulkfall.t.tree.data)
```

### Cyst Beta Diversity
```{r cystbeta, echo=FALSE}
pscyst.t.tree <- rtree(ntaxa(pscyst.t), rooted = TRUE, tip.label = taxa_names(pscyst.t))
#plot(pscyst.t.tree)
pscyst.t <- merge_phyloseq(pscyst.t, pscyst.t.tree)
```

#### Cyst Unweighted
```{r cystunweight, echo=FALSE}
pscysttunifrac_dist <- phyloseq::distance(pscyst.t, method="unifrac", weighted=F)
ordination = ordinate(pscyst.t, method="PCoA", distance=pscysttunifrac_dist)
#plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscyst.t.tree.data <- data.frame(sample_data(pscyst.t))
adonis2(pscysttunifrac_dist ~Plot*Season*Num_Cysts, pscyst.t.tree.data)
```

#### Cyst Weighted
```{r cystweight, echo=FALSE}
pscysttunifrac_dist <- phyloseq::distance(pscyst.t, method="unifrac", weighted=T)
ordination = ordinate(pscyst.t, method="PCoA", distance=pscysttunifrac_dist)
#plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscyst.t, ordination, color = "Plot", shape = "Season") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscyst.t.tree.data <- data.frame(sample_data(pscyst.t))
adonis2(pscysttunifrac_dist ~Plot*Season*Num_Cysts, pscyst.t.tree.data)
```

### Cyst Mid season Beta Diversity
```{r cystmidbeta, echo=FALSE}
pscystmid.t.tree <- rtree(ntaxa(pscystmid.t), rooted = TRUE, tip.label = taxa_names(pscystmid.t))
#plot(pscystmid.t.tree)
pscystmid.t <- merge_phyloseq(pscystmid.t, pscystmid.t.tree)
```

#### Unweighted Cysts Mid
```{r cystmidunweight, echo=FALSE}
pscystmidtunifrac_dist <- phyloseq::distance(pscystmid.t, method="unifrac", weighted=F)
ordination = ordinate(pscystmid.t, method="PCoA", distance=pscystmidtunifrac_dist)
#plot_ordination(pscystmid.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystmid.t.tree.data <- data.frame(sample_data(pscystmid.t))
adonis2(pscystmidtunifrac_dist ~Plot*Num_Cysts, pscystmid.t.tree.data)
```

#### Weighted Cysts Mid
```{r cystmidweight, echo=FALSE}
pscystmidtunifrac_dist <- phyloseq::distance(pscystmid.t, method="unifrac", weighted=T)
ordination = ordinate(pscystmid.t, method="PCoA", distance=pscystmidtunifrac_dist)
#plot_ordination(pscystmid.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystmid.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystmid.t.tree.data <- data.frame(sample_data(pscystmid.t))
adonis2(pscystmidtunifrac_dist ~Plot, pscystmid.t.tree.data)
```

### Cyst Fall season Beta Diversity
```{r cystfallbeta, echo=FALSE}
pscystfall.t.tree <- rtree(ntaxa(pscystfall.t), rooted = TRUE, tip.label = taxa_names(pscystfall.t))
#plot(pscystfall.t.tree)
pscystfall.t <- merge_phyloseq(pscystfall.t, pscystfall.t.tree)
```

#### Unweighted Cysts Fall
```{r cystfallunweight, echo=FALSE}
pscystfalltunifrac_dist <- phyloseq::distance(pscystfall.t, method="unifrac", weighted=F)
ordination = ordinate(pscystfall.t, method="PCoA", distance=pscystfalltunifrac_dist)
#plot_ordination(pscystfall.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystfall.t.tree.data <- data.frame(sample_data(pscystfall.t))
adonis2(pscystfalltunifrac_dist ~Plot, pscystfall.t.tree.data)
```

#### Weighted Cysts Fall
```{r cystfallweight, echo=FALSE}
pscystfalltunifrac_dist <- phyloseq::distance(pscystfall.t, method="unifrac", weighted=T)
ordination = ordinate(pscystfall.t, method="PCoA", distance=pscystfalltunifrac_dist)
#plot_ordination(pscystfall.t, ordination, color = "Plot", shape = "Season", label = "Num_Cysts") + theme(aspect.ratio=1) + geom_point(size=2)
plot_ordination(pscystfall.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
pscystfall.t.tree.data <- data.frame(sample_data(pscystfall.t))
adonis2(pscystfalltunifrac_dist ~Plot, pscystfall.t.tree.data)
```

### Root Beta Diversity
```{r rootbeta, echo=FALSE}
psroot.t.tree <- rtree(ntaxa(psroot.t), rooted = TRUE, tip.label = taxa_names(psroot.t))
#plot(psroot.t.tree)
psroot.t <- merge_phyloseq(psroot.t, psroot.t.tree)
```

#### Root Unweighted
```{r rootunweight, echo=FALSE}
psroottunifrac_dist <- phyloseq::distance(psroot.t, method="unifrac", weighted=F)
ordination = ordinate(psroot.t, method="PCoA", distance=psroottunifrac_dist)
plot_ordination(psroot.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psroot.t.tree.data <- data.frame(sample_data(psroot.t))
adonis2(psroottunifrac_dist ~Plot, psroot.t.tree.data)
```

#### Root Weighted
```{r rootweight, echo=FALSE}
psroottunifrac_dist <- phyloseq::distance(psroot.t, method="unifrac", weighted=T)
ordination = ordinate(psroot.t, method="PCoA", distance=psroottunifrac_dist)
plot_ordination(psroot.t, ordination, color = "Plot") + theme(aspect.ratio=1) + geom_point(size=2)
```

```{r}
psroot.t.tree.data <- data.frame(sample_data(psroot.t))
adonis2(psroottunifrac_dist ~Plot, psroot.t.tree.data)
```

## Bar plots highlighting most abundant taxa

### Rhizosphere Phyla present with at least 1% abundance
```{r rhizphyla, echo=FALSE}
psrhiz.rab <- transform_sample_counts(psrhiz, function(x) x / sum(x))
psrhiz.phylum <- tax_glom(psrhiz.rab, taxrank = "Phylum", NArm = FALSE)

rhiz.phylum <- psmelt(psrhiz.phylum)
rhiz.phylum$Phylum <- as.character(rhiz.phylum$Phylum)
max.rhiz.phylum <- ddply(rhiz.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.phylum <- max.rhiz.phylum[max.rhiz.phylum$max <= .01,]$Phylum
rhiz.phylum[rhiz.phylum$Phylum %in% maxremainder.rhiz.phylum,]$Phylum <- 'Other'

ggplot(rhiz.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Rhizosphere Classes present with at least 2% abundance
```{r rhizclass, echo=FALSE}
psrhiz.class <- tax_glom(psrhiz.rab, taxrank = "Class", NArm = FALSE)

rhiz.class <- psmelt(psrhiz.class)
rhiz.class$Class <- as.character(rhiz.class$Class)
max.rhiz.class <- ddply(rhiz.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.class <- max.rhiz.class[max.rhiz.class$max <= .01,]$Class
rhiz.class[rhiz.class$Class %in% maxremainder.rhiz.class,]$Class <- 'Other'

ggplot(rhiz.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Rhizosphere Families present with at least 3% abundance
```{r rhizfamily, echo=FALSE}
psrhiz.family <- tax_glom(psrhiz.rab, taxrank = "Family", NArm = FALSE)

rhiz.family <- psmelt(psrhiz.family)
rhiz.family$Family <- as.character(rhiz.family$Family)
max.rhiz.family <- ddply(rhiz.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.family <- max.rhiz.family[max.rhiz.family$max <= .02,]$Family
rhiz.family[rhiz.family$Family %in% maxremainder.rhiz.family,]$Family <- 'Other'

ggplot(rhiz.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Rhizosphere Genera present with at least 3% abundance
```{r rhizgenus, echo=FALSE}
psrhiz.genus <- tax_glom(psrhiz.rab, taxrank = "Genus", NArm = FALSE)

rhiz.genus <- psmelt(psrhiz.genus)
rhiz.genus$Genus <- as.character(rhiz.genus$Genus)
max.rhiz.genus <- ddply(rhiz.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.rhiz.genus <- max.rhiz.genus[max.rhiz.genus$max <= .03,]$Genus
rhiz.genus[rhiz.genus$Genus %in% maxremainder.rhiz.genus,]$Genus <- 'Other'

ggplot(rhiz.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot*Season, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()
) + guides(fill=guide_legend(title='Genus')   ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk Mid Phyla present with at least 1% abundance
```{r bulkmidphyla, echo=FALSE}
psbulkmid.rab <- transform_sample_counts(psbulkmid, function(x) x / sum(x))
psbulkmid.phylum <- tax_glom(psbulkmid.rab, taxrank = "Phylum", NArm = FALSE)

bulkmid.phylum <- psmelt(psbulkmid.phylum)
bulkmid.phylum$Phylum <- as.character(bulkmid.phylum$Phylum)
max.bulkmid.phylum <- ddply(bulkmid.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.phylum <- max.bulkmid.phylum[max.bulkmid.phylum$max <= .01,]$Phylum
bulkmid.phylum[bulkmid.phylum$Phylum %in% maxremainder.bulkmid.phylum,]$Phylum <- 'Other'

ggplot(bulkmid.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Bulk Mid Classes present with at least 2% abundance
```{r bulkmidclass, echo=FALSE}
psbulkmid.class <- tax_glom(psbulkmid.rab, taxrank = "Class", NArm = FALSE)

bulkmid.class <- psmelt(psbulkmid.class)
bulkmid.class$Class <- as.character(bulkmid.class$Class)
max.bulkmid.class <- ddply(bulkmid.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.class <- max.bulkmid.class[max.bulkmid.class$max <= .02,]$Class
bulkmid.class[bulkmid.class$Class %in% maxremainder.bulkmid.class,]$Class <- 'Other'

ggplot(bulkmid.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Class')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk Mid Families present with at least 5% abundance
```{r bulkmidfamily, echo=FALSE}
psbulkmid.family <- tax_glom(psbulkmid.rab, taxrank = "Family", NArm = FALSE)

bulkmid.family <- psmelt(psbulkmid.family)
bulkmid.family$Family <- as.character(bulkmid.family$Family)
max.bulkmid.family <- ddply(bulkmid.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.family <- max.bulkmid.family[max.bulkmid.family$max <= .05,]$Family
bulkmid.family[bulkmid.family$Family %in% maxremainder.bulkmid.family,]$Family <- 'Other'

ggplot(bulkmid.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Family')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

```

### Bulk Mid Genera present with at least 5% abundance
```{r bulkmidgenus, echo=FALSE}
psbulkmid.genus <- tax_glom(psbulkmid.rab, taxrank = "Genus", NArm = FALSE)

bulkmid.genus <- psmelt(psbulkmid.genus)
bulkmid.genus$Genus <- as.character(bulkmid.genus$Genus)
max.bulkmid.genus <- ddply(bulkmid.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.bulkmid.genus <- max.bulkmid.genus[max.bulkmid.genus$max <= .05,]$Genus
bulkmid.genus[bulkmid.genus$Genus %in% maxremainder.bulkmid.genus,]$Genus <- 'Other'

ggplot(bulkmid.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```



### Bulk Fall Phyla present with at least 1% abundance
```{r bulkfallphyla, echo=FALSE}
psbulkfall.rab <- transform_sample_counts(psbulkfall, function(x) x / sum(x))
psbulkfall.phylum <- tax_glom(psbulkfall.rab, taxrank = "Phylum", NArm = FALSE)

bulkfall.phylum <- psmelt(psbulkfall.phylum)
bulkfall.phylum$Phylum <- as.character(bulkfall.phylum$Phylum)
max.bulkfall.phylum <- ddply(bulkfall.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.phylum <- max.bulkfall.phylum[max.bulkfall.phylum$max <= .01,]$Phylum
bulkfall.phylum[bulkfall.phylum$Phylum %in% maxremainder.bulkfall.phylum,]$Phylum <- 'Other'

ggplot(bulkfall.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"),"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Bulk Fall Classes present with at least 2% abundance
```{r bulkfallclass, , echo=FALSE}
psbulkfall.class <- tax_glom(psbulkfall.rab, taxrank = "Class", NArm = FALSE)

bulkfall.class <- psmelt(psbulkfall.class)
bulkfall.class$Class <- as.character(bulkfall.class$Class)
max.bulkfall.class <- ddply(bulkfall.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.class <- max.bulkfall.class[max.bulkfall.class$max <= .02,]$Class
bulkfall.class[bulkfall.class$Class %in% maxremainder.bulkfall.class,]$Class <- 'Other'

ggplot(bulkfall.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk Fall families present with at least 5% abundance
```{r bulkfallfam, , echo=FALSE}
psbulkfall.family <- tax_glom(psbulkfall.rab, taxrank = "Family", NArm = FALSE)

bulkfall.family <- psmelt(psbulkfall.family)
bulkfall.family$Family <- as.character(bulkfall.family$Family)
max.bulkfall.family <- ddply(bulkfall.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.family <- max.bulkfall.family[max.bulkfall.family$max <= .05,]$Family
bulkfall.family[bulkfall.family$Family %in% maxremainder.bulkfall.family,]$Family <- 'Other'

ggplot(bulkfall.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk Fall Genera present with at least 5% abundance
```{r bulkfallgenera, , echo=FALSE}
psbulkfall.genus <- tax_glom(psbulkfall.rab, taxrank = "Genus", NArm = FALSE)

bulkfall.genus <- psmelt(psbulkfall.genus)
bulkfall.genus$Genus <- as.character(bulkfall.genus$Genus)
max.bulkfall.genus <- ddply(bulkfall.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.bulkfall.genus <- max.bulkfall.genus[max.bulkfall.genus$max <= .05,]$Genus
bulkfall.genus[bulkfall.genus$Genus %in% maxremainder.bulkfall.genus,]$Genus <- 'Other'

ggplot(bulkfall.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), 
"Other")))) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(strip.text = element_text(),axis.text.x = element_blank(), axis.ticks.x = element_blank(
), strip.background = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')  
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


### Cyst Mid season Phyla present with at least 1%
```{r cystmidphyla, echo=FALSE}
pscystmid.rab <- transform_sample_counts(pscystmid, function(x) x/sum(x))

pscystmid.phylum <- tax_glom(pscystmid.rab, taxrank = "Phylum", NArm = FALSE)
#plot_bar(pscystmid.phylum, x="Sample", fill="Phylum") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystmid.phylum <- psmelt(pscystmid.phylum)
# convert Phylum to a character vector from a factor because R
cystmid.phylum$Phylum <- as.character(cystmid.phylum$Phylum)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystmid.phylum <- ddply(cystmid.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.phylum <- max.cystmid.phylum[max.cystmid.phylum$max <= .01,]$Phylum
cystmid.phylum[cystmid.phylum$Phylum %in% maxremainder.cystmid.phylum,]$Phylum <- 'Other'
ggplot(cystmid.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Cyst Mid season Classes present with at least 1% abundance
```{r cystmidclass, echo=FALSE}
pscystmid.class <- tax_glom(pscystmid.rab, taxrank = "Class", NArm = FALSE)
#plot_bar(pscystmid.class, x="Sample", fill="Class") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystmid.class <- psmelt(pscystmid.class)
# convert Phylum to a character vector from a factor because R
cystmid.class$Class <- as.character(cystmid.class$Class)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystmid.class <- ddply(cystmid.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.class <- max.cystmid.class[max.cystmid.class$max <= .01,]$Class
cystmid.class[cystmid.class$Class %in% maxremainder.cystmid.class,]$Class <- 'Other'
ggplot(cystmid.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), 
          panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Mid season Families present with at least 5% abundance
```{r cystmidfam, echo=FALSE}
pscystmid.family <- tax_glom(pscystmid.rab, taxrank = "Family", NArm = FALSE)
#plot_bar(pscystmid.family, x="Sample", fill="Family") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystmid.family <- psmelt(pscystmid.family)
# convert taxa to a character vector from a factor because R
cystmid.family$Family <- as.character(cystmid.family$Family)
#Any taxa that is present in a sample greater than 1%
max.cystmid.family <- ddply(cystmid.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.family <- max.cystmid.family[max.cystmid.family$max <= .05,]$Family
cystmid.family[cystmid.family$Family %in% maxremainder.cystmid.family,]$Family <- 'Other'
ggplot(cystmid.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Mid season Genera present with at least 5% abundance
```{r cystmidgen, echo=FALSE}
pscystmid.genus <- tax_glom(pscystmid.rab, taxrank = "Genus", NArm = FALSE)
#plot_bar(pscystmid.genus, x="Sample", fill="Genus") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystmid.genus <- psmelt(pscystmid.genus)
# convert taxa to a character vector from a factor because R
cystmid.genus$Genus <- as.character(cystmid.genus$Genus)
#Any taxa that is present in a sample greater than 1%
max.cystmid.genus <- ddply(cystmid.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.cystmid.genus <- max.cystmid.genus[max.cystmid.genus$max <= .05,]$Genus
cystmid.genus[cystmid.genus$Genus %in% maxremainder.cystmid.genus,]$Genus <- 'Other'

ggplot(cystmid.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')
          ) +  scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Fall Phyla present with at least 1%
```{r cystfallphyla, echo=FALSE}
pscystfall.rab <- transform_sample_counts(pscystfall, function(x) x/sum(x))

pscystfall.phylum <- tax_glom(pscystfall.rab, taxrank = "Phylum", NArm = FALSE)
#plot_bar(pscystfall.phylum, x="Sample", fill="Phylum") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystfall.phylum <- psmelt(pscystfall.phylum)
# convert Phylum to a character vector from a factor because R
cystfall.phylum$Phylum <- as.character(cystfall.phylum$Phylum)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystfall.phylum <- ddply(cystfall.phylum, ~Phylum, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.phylum <- max.cystfall.phylum[max.cystfall.phylum$max <= .01,]$Phylum
cystfall.phylum[cystfall.phylum$Phylum %in% maxremainder.cystfall.phylum,]$Phylum <- 'Other'
ggplot(cystfall.phylum, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = c(setdiff(Phylum, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Phylum')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Cyst Fall Classes present with at least 1% abundance
```{r cystfallclass, echo=FALSE}
pscystfall.class <- tax_glom(pscystfall.rab, taxrank = "Class", NArm = FALSE)
#plot_bar(pscystfall.class, x="Sample", fill="Class") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystfall.class <- psmelt(pscystfall.class)
# convert Phylum to a character vector from a factor because R
cystfall.class$Class <- as.character(cystfall.class$Class)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.cystfall.class <- ddply(cystfall.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.class <- max.cystfall.class[max.cystfall.class$max <= .01,]$Class
cystfall.class[cystfall.class$Class %in% maxremainder.cystfall.class,]$Class <- 'Other'
ggplot(cystfall.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), 
          panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Fall Families present with at least 5% abundance
```{r cystfallfam, echo=FALSE}
pscystfall.family <- tax_glom(pscystfall.rab, taxrank = "Family", NArm = FALSE)
#plot_bar(pscystfall.family, x="Sample", fill="Family") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystfall.family <- psmelt(pscystfall.family)
# convert taxa to a character vector from a factor because R
cystfall.family$Family <- as.character(cystfall.family$Family)
#Any taxa that is present in a sample greater than 1%
max.cystfall.family <- ddply(cystfall.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.family <- max.cystfall.family[max.cystfall.family$max <= .05,]$Family
cystfall.family[cystfall.family$Family %in% maxremainder.cystfall.family,]$Family <- 'Other'
ggplot(cystfall.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cyst Fall Genera present with at least 5% abundance
```{r cystfallgen, echo=FALSE}
pscystfall.genus <- tax_glom(pscystfall.rab, taxrank = "Genus", NArm = FALSE)
#plot_bar(pscystfall.genus, x="Sample", fill="Genus") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank(), panel.background = element_blank())

cystfall.genus <- psmelt(pscystfall.genus)
# convert taxa to a character vector from a factor because R
cystfall.genus$Genus <- as.character(cystfall.genus$Genus)
#Any taxa that is present in a sample greater than 1%
max.cystfall.genus <- ddply(cystfall.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.cystfall.genus <- max.cystfall.genus[max.cystfall.genus$max <= .05,]$Genus
cystfall.genus[cystfall.genus$Genus %in% maxremainder.cystfall.genus,]$Genus <- 'Other'

ggplot(cystfall.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genus')
          ) +  scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Root Phyla
```{r rootphyla, echo=FALSE}
psroot.rab <- transform_sample_counts(psroot, function(x) x/sum(x))

psroot.phylum <- tax_glom(psroot.rab, taxrank = "Phylum", NArm = FALSE)
#plot_bar(psroot.phylum, x="Sample", fill="Phylum") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)  , strip.background = element_blank(), panel.background = element_blank())

root.phylum <- psmelt(psroot.phylum)
# convert Phylum to a character vector from a factor because R
root.phylum$Phylum <- as.character(root.phylum$Phylum)
ggplot(root.phylum, aes(x = Sample, y = Abundance, fill = Phylum))+ geom_bar(stat = "identity"
)+ facet_grid(~Plot, scales = "free", space = "free") + theme(axis.text.x = element_blank(), 
strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()
) + scale_fill_viridis(discrete = TRUE, option = "turbo")
```

### Root Classes present with at least 1% abundance
```{r rootclass, echo=FALSE}
psroot.class <- tax_glom(psroot.rab, taxrank = "Class", NArm = FALSE)
#plot_bar(psroot.class, x="Sample", fill="Class") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)  , strip.background = element_blank(), panel.background = element_blank())

#create dataframe from phyloseq object
root.class <- psmelt(psroot.class)
# convert Phylum to a character vector from a factor because R
root.class$Class <- as.character(root.class$Class)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.root.class <- ddply(root.class, ~Class, function(x) c(max=max(x$Abundance)))
maxremainder.root.class <- max.root.class[max.root.class$max <= .01,]$Class
root.class[root.class$Class %in% maxremainder.root.class,]$Class <- 'Other'

ggplot(root.class, aes(x = Sample, y = Abundance, fill = factor(Class, levels = c(setdiff(Class, "Other"), "Other")))
) + geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Class')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Root Families present with at least 3% abundance
```{r rootfam, echo=FALSE}
psroot.family <- tax_glom(psroot.rab, taxrank = "Family", NArm = FALSE)
#plot_bar(psroot.family, x="Sample", fill="Family") + facet_grid(~Plot, scales = "free", space = "free") + theme(strip.text = element_text(angle = 315),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)  , strip.background = element_blank(), panel.background = element_blank())

root.family <- psmelt(psroot.family)
# convert Phylum to a character vector from a factor because R
root.family$Family <- as.character(root.family$Family)
#Any taxa that is present in a sample greater than 1%
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.root.family <- ddply(root.family, ~Family, function(x) c(max=max(x$Abundance)))
maxremainder.root.family <- max.root.family[max.root.family$max <= .03,]$Family
root.family[root.family$Family %in% maxremainder.root.family,]$Family <- 'Other'

ggplot(root.family, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Family')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Root Genera present with at least 3% abundance
```{r rootgen, echo=FALSE}
psroot.genus <- tax_glom(psroot.rab, taxrank = "Genus", NArm=FALSE)

# create dataframe from phyloseq object
root.genus <- psmelt(psroot.genus)
# convert Phylum to a character vector from a factor because R
root.genus$Genus <- as.character(root.genus$Genus)
#must run the psmelt and as.character command before this. It will not update if the percentages are changed 
max.root.genus <- ddply(root.genus, ~Genus, function(x) c(max=max(x$Abundance)))
maxremainder.root.genus <- max.root.genus[max.root.genus$max <= .03,]$Genus
root.genus[root.genus$Genus %in% maxremainder.root.genus,]$Genus <- 'Other'

ggplot(root.genus, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity")+ facet_grid(~Plot, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
          ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

```


## Core Microbiome Plots

### All Samples Core Microbiome

#### Genera present with at least .01% abundance in 50% of all samples

```{r echo = FALSE}
psall <- ps3
psall.rab <- transform_sample_counts(psall, function(x) x/sum(x))
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .50*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### Genera present with at least .01% abundance in 60% of all samples

```{r echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .60*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


#### Genera present with at least .01% abundance in 70% of all samples

```{r echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .70*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```


#### Genera present with at least .01% abundance in 75% of all samples

```{r echo=FALSE}
filterall.rab <- phyloseq::genefilter_sample(psall.rab, filterfun_sample(function(x) x / sum(x) > .001), A = .75*nsamples(psall.rab))

allcore.rab <- prune_taxa(filterall.rab, psall.rab)
allcore.genus.rab <- tax_glom(allcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
allcore.genus.rab.df <- psmelt(allcore.genus.rab)
# convert Genus to a character vector from a factor because R
allcore.genus.rab.df$Genus <- as.character(allcore.genus.rab.df$Genus)
ggplot(allcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity", width = 1) + facet_grid(~SampleType*Season*SoilType, scales = "free", space = "free")+ theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank(), legend.text = element_text(size = 7)) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```



### Rhizosphere Core Microbiome

#### Subset rhizosphere samples by plot

```{r echo = FALSE}
psrhizSa <- subset_samples(psrhiz, Plot == "Sa")
psrhizS3 <- subset_samples(psrhiz, Plot == "S3")
psrhizSs <- subset_samples(psrhiz, Plot == "Ss")
psrhizSr <- subset_samples(psrhiz, Plot == "Sr")
```

#### SA Rhizosphere Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psrhizSa.rab <- transform_sample_counts(psrhizSa, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizSa.rab, filterfun_sample(function(x) x / sum(x) > .02))
rhizSacore.rab <- prune_taxa(filter.rab, psrhizSa.rab)
rhizSacore.genus.rab <- tax_glom(rhizSacore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizSacore.genus.rab.df <- psmelt(rhizSacore.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSacore.genus.rab.df$Genus <- as.character(rhizSacore.genus.rab.df$Genus)
ggplot(rhizSacore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Rhizosphere Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psrhizS3.rab <- transform_sample_counts(psrhizS3, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizS3.rab, filterfun_sample(function(x) x / sum(x) > .02))
rhizS3core.rab <- prune_taxa(filter.rab, psrhizS3.rab)
rhizS3core.genus.rab <- tax_glom(rhizS3core.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizS3core.genus.rab.df <- psmelt(rhizS3core.genus.rab)
# convert Genus to a character vector from a factor because R
rhizS3core.genus.rab.df$Genus <- as.character(rhizS3core.genus.rab.df$Genus)
ggplot(rhizS3core.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Rhizosphere Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psrhizSr.rab <- transform_sample_counts(psrhizSr, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizSr.rab, filterfun_sample(function(x) x / sum(x) > .02))
rhizSrcore.rab <- prune_taxa(filter.rab, psrhizSr.rab)
rhizSrcore.genus.rab <- tax_glom(rhizSrcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizSrcore.genus.rab.df <- psmelt(rhizSrcore.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSrcore.genus.rab.df$Genus <- as.character(rhizSrcore.genus.rab.df$Genus)
ggplot(rhizSrcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Rhizosphere Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psrhizSs.rab <- transform_sample_counts(psrhizSs, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrhizSs.rab, filterfun_sample(function(x) x / sum(x) > .02))
rhizSscore.rab <- prune_taxa(filter.rab, psrhizSs.rab)
rhizSscore.genus.rab <- tax_glom(rhizSscore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rhizSscore.genus.rab.df <- psmelt(rhizSscore.genus.rab)
# convert Genus to a character vector from a factor because R
rhizSscore.genus.rab.df$Genus <- as.character(rhizSscore.genus.rab.df$Genus)
ggplot(rhizSscore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
   panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
   ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```



### Bulk Soil Mid Season Core Microbiome

#### Subset bulk samples into plot and mid season

```{r echo = FALSE}
psbulkSa <- subset_samples(psbulk, Plot == "Sa")
psbulkS3 <- subset_samples(psbulk, Plot == "S3")
psbulkSs <- subset_samples(psbulk, Plot == "Ss")
psbulkSr <- subset_samples(psbulk, Plot == "Sr")
psbulkSamid <- subset_samples(psbulkSa, Season == "Mid")
psbulkS3mid <- subset_samples(psbulkS3, Season == "Mid")
psbulkSsmid <- subset_samples(psbulkSs, Season == "Mid")
psbulkSrmid <- subset_samples(psbulkSr, Season == "Mid")
```

#### SA Bulk Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psbulkSamid.rab <- transform_sample_counts(psbulkSamid, function(x) x/sum(x))
filterbSamid.rab <- phyloseq::genefilter_sample(psbulkSamid.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkSamidcore.rab <- prune_taxa(filterbSamid.rab, psbulkSamid.rab)
bulkSamidcore.genus.rab <- tax_glom(bulkSamidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSamidcore.genus.rab.df <- psmelt(bulkSamidcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSamidcore.genus.rab.df$Genus <- as.character(bulkSamidcore.genus.rab.df$Genus)
ggplot(bulkSamidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Bulk Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psbulkS3mid.rab <- transform_sample_counts(psbulkS3mid, function(x) x/sum(x))
filterbS3mid.rab <- phyloseq::genefilter_sample(psbulkS3mid.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkS3midcore.rab <- prune_taxa(filterbS3mid.rab, psbulkS3mid.rab)
bulkS3midcore.genus.rab <- tax_glom(bulkS3midcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkS3midcore.genus.rab.df <- psmelt(bulkS3midcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3midcore.genus.rab.df$Genus <- as.character(bulkS3midcore.genus.rab.df$Genus)
ggplot(bulkS3midcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Bulk Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psbulkSrmid.rab <- transform_sample_counts(psbulkSrmid, function(x) x/sum(x))
filterbSrmid.rab <- phyloseq::genefilter_sample(psbulkSrmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkSrmidcore.rab <- prune_taxa(filterbSrmid.rab, psbulkSrmid.rab)
bulkSrmidcore.genus.rab <- tax_glom(bulkSrmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSrmidcore.genus.rab.df <- psmelt(bulkSrmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSrmidcore.genus.rab.df$Genus <- as.character(bulkSrmidcore.genus.rab.df$Genus)
ggplot(bulkSrmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Bulk Soil Mid season core genera (2% abundance or higher)

```{r echo = FALSE}
psbulkSsmid.rab <- transform_sample_counts(psbulkSsmid, function(x) x/sum(x))
filterbSsmid.rab <- phyloseq::genefilter_sample(psbulkSsmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkSsmidcore.rab <- prune_taxa(filterbSsmid.rab, psbulkSsmid.rab)
bulkSsmidcore.genus.rab <- tax_glom(bulkSsmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSsmidcore.genus.rab.df <- psmelt(bulkSsmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSsmidcore.genus.rab.df$Genus <- as.character(bulkSsmidcore.genus.rab.df$Genus)
ggplot(bulkSsmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Bulk Soil Fall Season Core Microbiome

#### Subset bulk samples into fall season

```{r echo = FALSE}
psbulkSafall <- subset_samples(psbulkSa, Season == "Fall")
psbulkS3fall <- subset_samples(psbulkS3, Season == "Fall")
psbulkSsfall <- subset_samples(psbulkSs, Season == "Fall")
psbulkSrfall <- subset_samples(psbulkSr, Season == "Fall")
```

#### SA Bulk Soil Fall season core genera (2% abundance or higher)

```{r echo = FALSE}
psbulkSafall.rab <- transform_sample_counts(psbulkSafall, function(x) x/sum(x))
filterbSafall.rab <- phyloseq::genefilter_sample(psbulkSafall.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkSafallcore.rab <- prune_taxa(filterbSafall.rab, psbulkSafall.rab)
bulkSafallcore.genus.rab <- tax_glom(bulkSafallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSafallcore.genus.rab.df <- psmelt(bulkSafallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSafallcore.genus.rab.df$Genus <- as.character(bulkSafallcore.genus.rab.df$Genus)
ggplot(bulkSafallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Bulk Soil Fall season core genera (2% abundance or higher)

```{r  echo=FALSE}
psbulkS3fall.rab <- transform_sample_counts(psbulkS3fall, function(x) x/sum(x))
filterbS3fall.rab <- phyloseq::genefilter_sample(psbulkS3fall.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkS3fallcore.rab <- prune_taxa(filterbS3fall.rab, psbulkS3fall.rab)
bulkS3fallcore.genus.rab <- tax_glom(bulkS3fallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkS3fallcore.genus.rab.df <- psmelt(bulkS3fallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkS3fallcore.genus.rab.df$Genus <- as.character(bulkS3fallcore.genus.rab.df$Genus)
ggplot(bulkS3fallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Bulk Soil Fall season core genera (2% abundance or higher)

```{r echo=FALSE}
psbulkSrfall.rab <- transform_sample_counts(psbulkSrfall, function(x) x/sum(x))
filterbSrfall.rab <- phyloseq::genefilter_sample(psbulkSrfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkSrfallcore.rab <- prune_taxa(filterbSrfall.rab, psbulkSrfall.rab)
bulkSrfallcore.genus.rab <- tax_glom(bulkSrfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
bulkSrfallcore.genus.rab.df <- psmelt(bulkSrfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
bulkSrfallcore.genus.rab.df$Genus <- as.character(bulkSrfallcore.genus.rab.df$Genus)
ggplot(bulkSrfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Bulk Soil Fall season core genera (2% abundance or higher)

```{r echo=FALSE}
psbulkSsfall.rab <- transform_sample_counts(psbulkSsfall, function(x) x/sum(x))
filterbSsfall.rab <- phyloseq::genefilter_sample(psbulkSsfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
bulkSsfallcore.rab <- prune_taxa(filterbSsfall.rab, psbulkSsfall.rab)
bulkSsfallcore.genus.rab <- tax_glom(bulkSsfallcore.rab, taxrank = "Genus", NArm = FALSE)

bulkSsfallcore.genus.rab.df <- psmelt(bulkSsfallcore.genus.rab)

bulkSsfallcore.genus.rab.df$Genus <- as.character(bulkSsfallcore.genus.rab.df$Genus)
ggplot(bulkSsfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts Mid Season Core Microbiome

#### Subset cyst samples into plots & mid season

```{r echo=FALSE}
pscystSa <- subset_samples(pscyst, Plot == "Sa")
pscystS3 <- subset_samples(pscyst, Plot == "S3")
pscystSs <- subset_samples(pscyst, Plot == "Ss")
pscystSr <- subset_samples(pscyst, Plot == "Sr")
pscystSamid <- subset_samples(pscystSa, Season == "Mid")
pscystS3mid <- subset_samples(pscystS3, Season == "Mid")
pscystSsmid <- subset_samples(pscystSs, Season == "Mid")
pscystSrmid <- subset_samples(pscystSr, Season == "Mid")
```

#### SA Cyst Mid season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystSamid.rab <- transform_sample_counts(pscystSamid, function(x) x/sum(x))
filtercSamid.rab <- phyloseq::genefilter_sample(pscystSamid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSamidcore.rab <- prune_taxa(filtercSamid.rab, pscystSamid.rab)
cystSamidcore.genus.rab <- tax_glom(cystSamidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSamidcore.genus.rab.df <- psmelt(cystSamidcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSamidcore.genus.rab.df$Genus <- as.character(cystSamidcore.genus.rab.df$Genus)
ggplot(cystSamidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Cyst Mid season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystS3mid.rab <- transform_sample_counts(pscystS3mid, function(x) x/sum(x))
filtercS3mid.rab <- phyloseq::genefilter_sample(pscystS3mid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3midcore.rab <- prune_taxa(filtercS3mid.rab, pscystS3mid.rab)
cystS3midcore.genus.rab <- tax_glom(cystS3midcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystS3midcore.genus.rab.df <- psmelt(cystS3midcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3midcore.genus.rab.df$Genus <- as.character(cystS3midcore.genus.rab.df$Genus)
ggplot(cystS3midcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Cyst Mid season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystSrmid.rab <- transform_sample_counts(pscystSrmid, function(x) x/sum(x))
filtercSrmid.rab <- phyloseq::genefilter_sample(pscystSrmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSrmidcore.rab <- prune_taxa(filtercSrmid.rab, pscystSrmid.rab)
cystSrmidcore.genus.rab <- tax_glom(cystSrmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSrmidcore.genus.rab.df <- psmelt(cystSrmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSrmidcore.genus.rab.df$Genus <- as.character(cystSrmidcore.genus.rab.df$Genus)
ggplot(cystSrmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Cyst Mid season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystSsmid.rab <- transform_sample_counts(pscystSsmid, function(x) x/sum(x))
filtercSsmid.rab <- phyloseq::genefilter_sample(pscystSsmid.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSsmidcore.rab <- prune_taxa(filtercSsmid.rab, pscystSsmid.rab)
cystSsmidcore.genus.rab <- tax_glom(cystSsmidcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSsmidcore.genus.rab.df <- psmelt(cystSsmidcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSsmidcore.genus.rab.df$Genus <- as.character(cystSsmidcore.genus.rab.df$Genus)
ggplot(cystSsmidcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts Mid Season Core- Family

#### SA Cyst Mid season core family (2% abundance or higher)

```{r echo=FALSE}
cystSamidcore.Family.rab <- tax_glom(cystSamidcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystSamidcore.Family.rab.df <- psmelt(cystSamidcore.Family.rab)
# convert Family to a character vector from a factor because R
cystSamidcore.Family.rab.df$Family <- as.character(cystSamidcore.Family.rab.df$Family)
ggplot(cystSamidcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Cyst Mid season core family (2% abundance or higher)

```{r echo=FALSE}
cystS3midcore.Family.rab <- tax_glom(cystS3midcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystS3midcore.Family.rab.df <- psmelt(cystS3midcore.Family.rab)
# convert Family to a character vector from a factor because R
cystS3midcore.Family.rab.df$Family <- as.character(cystS3midcore.Family.rab.df$Family)
ggplot(cystS3midcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Cyst Mid season core family (2% abundance or higher)

```{r echo=FALSE}
cystSrmidcore.Family.rab <- tax_glom(cystSrmidcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystSrmidcore.Family.rab.df <- psmelt(cystSrmidcore.Family.rab)
# convert Family to a character vector from a factor because R
cystSrmidcore.Family.rab.df$Family <- as.character(cystSrmidcore.Family.rab.df$Family)
ggplot(cystSrmidcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Cyst Mid season core family (2% abundance or higher)

```{r echo=FALSE}
cystSsmidcore.Family.rab <- tax_glom(cystSsmidcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystSsmidcore.Family.rab.df <- psmelt(cystSsmidcore.Family.rab)
# convert Family to a character vector from a factor because R
cystSsmidcore.Family.rab.df$Family <- as.character(cystSsmidcore.Family.rab.df$Family)
ggplot(cystSsmidcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts Fall Season Core Microbiome

#### Subset cyst samples into fall season

```{r echo=FALSE}
pscystSafall <- subset_samples(pscystSa, Season == "Fall")
pscystS3fall <- subset_samples(pscystS3, Season == "Fall")
pscystSsfall <- subset_samples(pscystSs, Season == "Fall")
pscystSrfall <- subset_samples(pscystSr, Season == "Fall")
```

#### SA Cyst Fall season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystSafall.rab <- transform_sample_counts(pscystSafall, function(x) x/sum(x))
filtercSafall.rab <- phyloseq::genefilter_sample(pscystSafall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSafallcore.rab <- prune_taxa(filtercSafall.rab, pscystSafall.rab)
cystSafallcore.genus.rab <- tax_glom(cystSafallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSafallcore.genus.rab.df <- psmelt(cystSafallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSafallcore.genus.rab.df$Genus <- as.character(cystSafallcore.genus.rab.df$Genus)
ggplot(cystSafallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Cyst Fall season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystS3fall.rab <- transform_sample_counts(pscystS3fall, function(x) x/sum(x))
filtercS3fall.rab <- phyloseq::genefilter_sample(pscystS3fall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystS3fallcore.rab <- prune_taxa(filtercS3fall.rab, pscystS3fall.rab)
cystS3fallcore.genus.rab <- tax_glom(cystS3fallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystS3fallcore.genus.rab.df <- psmelt(cystS3fallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystS3fallcore.genus.rab.df$Genus <- as.character(cystS3fallcore.genus.rab.df$Genus)
ggplot(cystS3fallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Cyst Fall season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystSrfall.rab <- transform_sample_counts(pscystSrfall, function(x) x/sum(x))
filtercSrfall.rab <- phyloseq::genefilter_sample(pscystSrfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSrfallcore.rab <- prune_taxa(filtercSrfall.rab, pscystSrfall.rab)
cystSrfallcore.genus.rab <- tax_glom(cystSrfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSrfallcore.genus.rab.df <- psmelt(cystSrfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSrfallcore.genus.rab.df$Genus <- as.character(cystSrfallcore.genus.rab.df$Genus)
ggplot(cystSrfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Cyst Fall season core genera (2% abundance or higher)

```{r echo=FALSE}
pscystSsfall.rab <- transform_sample_counts(pscystSsfall, function(x) x/sum(x))
filtercSsfall.rab <- phyloseq::genefilter_sample(pscystSsfall.rab, filterfun_sample(function(x) x / sum(x) > .02))
cystSsfallcore.rab <- prune_taxa(filtercSsfall.rab, pscystSsfall.rab)
cystSsfallcore.genus.rab <- tax_glom(cystSsfallcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
cystSsfallcore.genus.rab.df <- psmelt(cystSsfallcore.genus.rab)
# convert Genus to a character vector from a factor because R
cystSsfallcore.genus.rab.df$Genus <- as.character(cystSsfallcore.genus.rab.df$Genus)
ggplot(cystSsfallcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Cysts Fall Season Core- Family

#### SA Cyst Fall season core family (2% abundance or higher)

```{r echo=FALSE}
cystSafallcore.Family.rab <- tax_glom(cystSafallcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystSafallcore.Family.rab.df <- psmelt(cystSafallcore.Family.rab)
# convert Family to a character vector from a factor because R
cystSafallcore.Family.rab.df$Family <- as.character(cystSafallcore.Family.rab.df$Family)
ggplot(cystSafallcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")

```

#### S3 Cyst Fall season core family (2% abundance or higher)

```{r echo=FALSE}
cystS3fallcore.Family.rab <- tax_glom(cystS3fallcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystS3fallcore.Family.rab.df <- psmelt(cystS3fallcore.Family.rab)
# convert Family to a character vector from a factor because R
cystS3fallcore.Family.rab.df$Family <- as.character(cystS3fallcore.Family.rab.df$Family)
ggplot(cystS3fallcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Cyst Fall season core family (2% abundance or higher)

```{r echo=FALSE}
cystSrfallcore.Family.rab <- tax_glom(cystSrfallcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystSrfallcore.Family.rab.df <- psmelt(cystSrfallcore.Family.rab)
# convert Family to a character vector from a factor because R
cystSrfallcore.Family.rab.df$Family <- as.character(cystSrfallcore.Family.rab.df$Family)
ggplot(cystSrfallcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Cyst Fall season core family (2% abundance or higher)

```{r echo=FALSE}
cystSsfallcore.Family.rab <- tax_glom(cystSsfallcore.rab, taxrank = "Family", NArm = FALSE)
#make dataframe
cystSsfallcore.Family.rab.df <- psmelt(cystSsfallcore.Family.rab)
# convert Family to a character vector from a factor because R
cystSsfallcore.Family.rab.df$Family <- as.character(cystSsfallcore.Family.rab.df$Family)
ggplot(cystSsfallcore.Family.rab.df, aes(x = Sample, y = Abundance, fill = factor(Family, levels = c(setdiff(Family, "Other"), "Other")))
) +  geom_bar(stat = "identity") + facet_grid(~Season, scales = "free", space = "free"
) + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(), panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

### Root core genera microbiome

#### Split root samples by plot

```{r echo=FALSE}
psrootSa <- subset_samples(psroot, Plot == "Sa")
psrootS3 <- subset_samples(psroot, Plot == "S3")
psrootSs <- subset_samples(psroot, Plot == "Ss")
psrootSr <- subset_samples(psroot, Plot == "Sr")
```

#### SA Root core genera (1% abundance or higher)

```{r echo=FALSE}
psrootSa.rab <- transform_sample_counts(psrootSa, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootSa.rab, filterfun_sample(function(x) x / sum(x) > .01))
rootSacore.rab <- prune_taxa(filter.rab, psrootSa.rab)
rootSacore.genus.rab <- tax_glom(rootSacore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootSacore.genus.rab.df <- psmelt(rootSacore.genus.rab)
# convert Genus to a character vector from a factor because R
rootSacore.genus.rab.df$Genus <- as.character(rootSacore.genus.rab.df$Genus)
ggplot(rootSacore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### S3 Root core genera (1% abundance or higher)

```{r echo=FALSE}
psrootS3.rab <- transform_sample_counts(psrootS3, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootS3.rab, filterfun_sample(function(x) x / sum(x) > .01))
rootS3core.rab <- prune_taxa(filter.rab, psrootS3.rab)
rootS3core.genus.rab <- tax_glom(rootS3core.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootS3core.genus.rab.df <- psmelt(rootS3core.genus.rab)
# convert Genus to a character vector from a factor because R
rootS3core.genus.rab.df$Genus <- as.character(rootS3core.genus.rab.df$Genus)
ggplot(rootS3core.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SR Root core genera (1% abundance or higher)

```{r echo=FALSE}
psrootSr.rab <- transform_sample_counts(psrootSr, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootSr.rab, filterfun_sample(function(x) x / sum(x) > .01))
rootSrcore.rab <- prune_taxa(filter.rab, psrootSr.rab)
rootSrcore.genus.rab <- tax_glom(rootSrcore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootSrcore.genus.rab.df <- psmelt(rootSrcore.genus.rab)
# convert Genus to a character vector from a factor because R
rootSrcore.genus.rab.df$Genus <- as.character(rootSrcore.genus.rab.df$Genus)
ggplot(rootSrcore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```

#### SS Root core genera (1% abundance or higher)

```{r echo=FALSE}
psrootSs.rab <- transform_sample_counts(psrootSs, function(x) x/sum(x))
filter.rab <- phyloseq::genefilter_sample(psrootSs.rab, filterfun_sample(function(x) x / sum(x) > .01))
rootSscore.rab <- prune_taxa(filter.rab, psrootSs.rab)
rootSscore.genus.rab <- tax_glom(rootSscore.rab, taxrank = "Genus", NArm = FALSE)
#make dataframe
rootSscore.genus.rab.df <- psmelt(rootSscore.genus.rab)
# convert Genus to a character vector from a factor because R
rootSscore.genus.rab.df$Genus <- as.character(rootSscore.genus.rab.df$Genus)
ggplot(rootSscore.genus.rab.df, aes(x = Sample, y = Abundance, fill = factor(Genus, levels = c(setdiff(Genus, "Other"), "Other")))) +  geom_bar(stat = "identity") + theme(axis.text.x = element_blank(), strip.background = element_blank(), axis.ticks.x = element_blank(),
  panel.background = element_blank()) + guides(fill=guide_legend(title='Genera')
  ) + scale_fill_viridis(discrete = TRUE, option = "turbo", na.value = "#A9A9A9")
```
