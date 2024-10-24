## Emily A Green, PhD

## Code for NetCoMi correlation networks that include both 16S and ITS data. 

### Exports of the 16S ps4 otu and tax table and the ITS ps3 otu and tax table were used.
#16Sps4_otutable.csv & 16Sps4_taxtable.csv

#ITSps3_otutable.csv & ITSps3_taxtable.csv

### In LTR16S_Metagenome.RData
#write.csv(otu_table(ps4), "16Sps4_otutable.csv")

#write.csv(tax_table(ps4), "16Sps4_taxtable.csv")

### In LTR_ITS_S3ARSWorkspace.RData
#write.csv(otu_table(ps3), "ITSps3_otutable.csv")

#write.csv(tax_table(ps3), "ITSps3_taxtable.csv")

### In excel, I merged the Merging the 16S and ITS data together to make new OTU & Tax tables, and metadata. 
#Samples were renamed in the metadatafile & Several samples were removed because they only appeared in one dataset

#remove EAG000054xITS

#remove EAG000084x16S

#remove EAG00089x16S

#remove EAG00199x16SFallCysts

#remove EAG00201x16SFallCysts

#Remove EAG00269xITSHornets

## Load packages 
library(dada2)

library(phyloseq)

library(ggplot2)

library(plyr)

library(vegan)

theme_set(theme_bw())

library(NetCoMi)

## Read in OTU table 
otu_table <- read.csv("Merged16SITS_otutable.csv")

otu_table.df <- data.frame(otu_table, row.names = 1)

otu_table.df <- otu_table(otu_table.df, taxa_are_rows = FALSE)

## Read in Tax table 
tax_table.df <- read.csv("Merged16SITS_taxtable.csv")

tax_table.df <- data.frame(tax_table.df, row.names = 1)

tax_table.matrix <- as.matrix(tax_table.df)

tax_table.matrix <- tax_table(tax_table.matrix)

## Read in metadata 
meta <- read.csv("Merged16SITSMetadata.csv", header = TRUE, row.names = 1)

## Make phyloseq object 
ps <- phyloseq(otu_table(otu_table.df, taxa_are_rows = FALSE), sample_data(meta), tax_table(tax_table.matrix))

### Split phyloseq into sample types 
psbulk <- subset_samples(ps, Sample_Type == "Bulk soil")

psrhiz <- subset_samples(ps, Sample_Type == "Rhizosphere")

pscyst <- subset_samples(ps, Sample_Type == "Cyst")

psroot <- subset_samples(ps, Sample_Type == "Roots")

### Remove taxa/OTU with zero reads 
psbulk <- prune_taxa(taxa_sums(psbulk) > 0, psbulk)

psrhiz <- prune_taxa(taxa_sums(psrhiz) > 0, psrhiz)

pscyst <- prune_taxa(taxa_sums(pscyst) > 0, pscyst)

psroot <- prune_taxa(taxa_sums(psroot) > 0, psroot)

### Transform abundance 
psbulk.rab <- transform_sample_counts(psbulk, function(x) x/sum(x))

psrhiz.rab <- transform_sample_counts(psrhiz, function(x) x/sum(x))

pscyst.rab <- transform_sample_counts(pscyst, function(x) x/sum(x))

psroot.rab <- transform_sample_counts(psroot, function(x) x/sum(x))

## NetCoMi- cysts 
filterpscyst.rab <- phyloseq::genefilter_sample(pscyst.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

pscyst.rab.filt <- prune_taxa(filterpscyst.rab, pscyst.rab)

pscyst.rab.glom <- tax_glom(pscyst.rab.filt, taxrank = "Genus", NArm = FALSE)

pscyst.rab.rename <-  renameTaxa(pscyst.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                     numDupli = "Genus", numDupliPat ="<name><.><num>")

net_pscystglom_spiec <- netConstruct(pscyst.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_pscystglom_spiec <- netAnalyze(net_pscystglom_spiec)

summary(analyze_pscystglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_pscystglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS Cysts glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as Cyst 16S ITS NetCoMi 1percent.pdf


## NetCoMi bulk soil 
filterpsbulk.rab <- phyloseq::genefilter_sample(psbulk.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

psbulk.rab.filt <- prune_taxa(filterpsbulk.rab, psbulk.rab)

psbulk.rab.glom <- tax_glom(psbulk.rab.filt, taxrank = "Genus", NArm = FALSE)

psbulk.rab.rename <-  renameTaxa(psbulk.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                 numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psbulkglom_spiec <- netConstruct(psbulk.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psbulkglom_spiec <- netAnalyze(net_psbulkglom_spiec)

summary(analyze_psbulkglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psbulkglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS Bulk glom > 1% abund"
 , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "Bulk 16S ITS NetCoMi 1percent.pdf"


## NetCoMi rhizosphere 
filterpsrhiz.rab <- phyloseq::genefilter_sample(psrhiz.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

psrhiz.rab.filt <- prune_taxa(filterpsrhiz.rab, psrhiz.rab)

psrhiz.rab.glom <- tax_glom(psrhiz.rab.filt, taxrank = "Genus", NArm = FALSE)

psrhiz.rab.rename <-  renameTaxa(psrhiz.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                 numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psrhizglom_spiec <- netConstruct(psrhiz.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psrhizglom_spiec <- netAnalyze(net_psrhizglom_spiec)

summary(analyze_psrhizglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psrhizglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS Rhiz glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "Rhiz 16S ITS NetCoMi 1percent.pdf"


## NetCoMi roots  
filterpsroot.rab <- phyloseq::genefilter_sample(psroot.rab, filterfun_sample(function(x) x / sum(x) > 0.02))

psroot.rab.filt <- prune_taxa(filterpsroot.rab, psroot.rab)

psroot.rab.glom <- tax_glom(psroot.rab.filt, taxrank = "Genus", NArm = FALSE)

psroot.rab.rename <-  renameTaxa(psroot.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                                 numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psrootglom_spiec <- netConstruct(psroot.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                     measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                     sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psrootglom_spiec <- netAnalyze(net_psrootglom_spiec)

summary(analyze_psrootglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psrootglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS Root glom > 2% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "Root 16S ITS NetCoMi 2percent.pdf"


## Split into treatments, but include all sample types 
psS3 <- subset_samples(ps, Plot == "S3")

psSa <- subset_samples(ps, Plot == "Sa")

psSr <- subset_samples(ps, Plot == "Sr")

psSs <- subset_samples(ps, Plot == "Ss")


### Remove taxa/OTU with zero reads 
psS3 <- prune_taxa(taxa_sums(psS3) > 0, psS3)

psSa <- prune_taxa(taxa_sums(psSa) > 0, psSa)

psSr <- prune_taxa(taxa_sums(psSr) > 0, psSr)

psSs <- prune_taxa(taxa_sums(psSs) > 0, psSs)


### Transform abundance 
psS3.rab <- transform_sample_counts(psS3, function(x) x/sum(x))

psSa.rab <- transform_sample_counts(psSa, function(x) x/sum(x))

psSr.rab <- transform_sample_counts(psSr, function(x) x/sum(x))

psSs.rab <- transform_sample_counts(psSs, function(x) x/sum(x))


## NetCoMi S3 
filterpsS3.rab <- phyloseq::genefilter_sample(psS3.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

psS3.rab.filt <- prune_taxa(filterpsS3.rab, psS3.rab)

psS3.rab.glom <- tax_glom(psS3.rab.filt, taxrank = "Genus", NArm = FALSE)

psS3.rab.rename <-  renameTaxa(psS3.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                               numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psS3glom_spiec <- netConstruct(psS3.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                   measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                   sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psS3glom_spiec <- netAnalyze(net_psS3glom_spiec)

summary(analyze_psS3glom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psS3glom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS S3 glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "S3 16S ITS NetCoMi 1percent.pdf" 


## NetCoMi Sa 
filterpsSa.rab <- phyloseq::genefilter_sample(psSa.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

psSa.rab.filt <- prune_taxa(filterpsSa.rab, psSa.rab)

psSa.rab.glom <- tax_glom(psSa.rab.filt, taxrank = "Genus", NArm = FALSE)

psSa.rab.rename <-  renameTaxa(psSa.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                               numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psSaglom_spiec <- netConstruct(psSa.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                   measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                   sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psSaglom_spiec <- netAnalyze(net_psSaglom_spiec)

summary(analyze_psSaglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psSaglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS SA glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "SA 16S ITS NetCoMi 1percent.pdf" 


## NetCoMi Sr 
filterpsSr.rab <- phyloseq::genefilter_sample(psSr.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

psSr.rab.filt <- prune_taxa(filterpsSr.rab, psSr.rab)

psSr.rab.glom <- tax_glom(psSr.rab.filt, taxrank = "Genus", NArm = FALSE)

psSr.rab.rename <-  renameTaxa(psSr.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                               numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psSrglom_spiec <- netConstruct(psSr.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                   measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                   sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psSrglom_spiec <- netAnalyze(net_psSrglom_spiec)

summary(analyze_psSrglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psSrglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS Sr glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "SR 16S ITS NetCoMi 1percent.pdf" 


## NetCoMi Ss 
filterpsSs.rab <- phyloseq::genefilter_sample(psSs.rab, filterfun_sample(function(x) x / sum(x) > 0.01))

psSs.rab.filt <- prune_taxa(filterpsSs.rab, psSs.rab)

psSs.rab.glom <- tax_glom(psSs.rab.filt, taxrank = "Genus", NArm = FALSE)

psSs.rab.rename <-  renameTaxa(psSs.rab.glom, pat = "<name>", substPat = "<name>_<subst_name>(<subst_R>)",
                               numDupli = "Genus", numDupliPat ="<name><.><num>")

net_psSsglom_spiec <- netConstruct(psSs.rab.rename, measure = "spieceasi", taxRank = "Genus", 
                                   measurePar = list(method = "mb", pulsar.params = list(rep.num = 50)),
                                   sparsMethod = "none", normMethod = "none", verbose = 3)

analyze_psSsglom_spiec <- netAnalyze(net_psSsglom_spiec)

summary(analyze_psSsglom_spiec)

#bold circles are hub nodes hub nodes = eigenvector centrality value above the empirical 95% quantile of all eigenvector centralities in the network 

plot(analyze_psSsglom_spiec, labelScale = TRUE, cexNodes = 2, labelFont = 3, rmSingles = TRUE, title1 = " 16S/ITS SS glom > 1% abund"
     , showTitle = TRUE, highlightHubs = TRUE, hubBorderWidth = 3)

#Saved as "SS 16S ITS NetCoMi 1percent.pdf" 
