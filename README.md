# SIFIC
# Immunology of Severe Febrile Illness in Children in the COVID-19 Era


Gene-set enrichment analysis was undertaken on previously published data from Jackson et al. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10312302/) that are available at ArrayExpress under accession E-MTAB-11671 and E-MTAB-12793. We used the R (version 4.4) implementation of gene-set enrichment analysis (fgsea, version 1.3) and the C7 immunologic signature gene sets available from as part of the Human MSigDB Collections (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp). We used Gene Ontology (GO) enrichment analysis. We compared gene enrichment between children with MIS-C and healthy controls, with definite bacterial infection and healthy controls, definite viral infection and healthy controls and children with MIS-C and definite bacterial infection. Figures were drawn in ggplot2.

## https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

library(fgsea)
library(dplyr)
library(goseq)
library(ggplot2)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

load("./DESeq_object_MISC_signature.RData")
load("./human_c7_v5p2.rdata")

design(dds_samples) <- ~ phenotype

#### MIS-C versus HC ####
dds_samples_new <- DESeq(dds_samples)
res_misc_hc <- results(dds_samples_new, contrast = c("phenotype", "HC", "MIS_C"))
res_misc_hc

res_misc_hc$GeneID <- row.names(res_misc_hc)

# Create ranks
gseaDat_misc_hc <- res_misc_hc
gseaDat_misc_hc$Entrez = mapIds(org.Hs.eg.db,
                                keys=row.names(gseaDat_misc_hc),
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first") #Ensembl gene IDs to Entrez gene IDs

gseaDat_misc_hc <- gseaDat_misc_hc[!is.na(gseaDat_misc_hc$Entrez),]

gseaDat_misc_hc$Symbol = mapIds(org.Hs.eg.db,
                                keys=gseaDat_misc_hc$Entrez,
                                column="SYMBOL",
                                keytype="ENTREZID",
                                multiVals="first") #Entrez gene IDs to gene symbols

ranks_misc_hc <- gseaDat_misc_hc$log2FoldChange
names(ranks_misc_hc) <- gseaDat_misc_hc$Entrez
head(ranks_misc_hc)
tail(ranks_misc_hc)

# Plot ranks
barplot(sort(ranks_misc_hc, decreasing = T))

# Load pathways
pathwaysH <- Hs.c7

# Conduct analysis
fgseaRes_misc_hc <- fgsea(pathwaysH, ranks_misc_hc, minSize = 15, maxSize = 500)

# Examining the top 10 results
head(fgseaRes_misc_hc)

# GSEA table plot
topUp <- fgseaRes_misc_hc %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes_misc_hc %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway],
              ranks_misc_hc, 
              fgseaRes_misc_hc, 
              gseaParam = 0.1)

#supportedOrganisms() %>% filter(str_detect(Genome, "hs"))

# List of differentially expressed genes
isSigGene <- gseaDat_misc_hc$padj < 0.01 & !is.na(gseaDat_misc_hc$padj)
genes <- as.integer(isSigGene)
names(genes) <- gseaDat_misc_hc$GeneID

# Probability weighting function
pwf <- nullp(genes, "hg19", "ensGene", bias.data = gseaDat_misc_hc$medianTxLength)

# GO enrichment analysis
goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:BP"))

# Plot the top 10
g.go.misc_hc <- goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "MIS-C versus pHC") + theme_bw()

ggsave("g.go.misc_hc.pdf", g.go.misc_hc, dpi = 300, height = 6, width = 8)

# KEGG enrichment analysis

#### DB versus HC ####
# dds_samples_new_db_hc <- dds_samples
# dds_samples_new_db_hc <- DESeq(dds_samples_new_db_hc)
res_db_hc <- results(dds_samples_new, contrast = c("phenotype", "HC", "DB"))
res_db_hc

res_db_hc$GeneID <- row.names(res_db_hc)

# Create ranks
gseaDat_db_hc <- res_db_hc
gseaDat_db_hc$Entrez = mapIds(org.Hs.eg.db,
                              keys=row.names(gseaDat_db_hc),
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first") #Ensembl gene IDs to Entrez gene IDs

gseaDat_db_hc <- gseaDat_db_hc[!is.na(gseaDat_db_hc$Entrez),]

gseaDat_db_hc$Symbol = mapIds(org.Hs.eg.db,
                              keys=gseaDat_db_hc$Entrez,
                              column="SYMBOL",
                              keytype="ENTREZID",
                              multiVals="first") #Entrez gene IDs to gene symbols

ranks_db_hc <- gseaDat_db_hc$log2FoldChange
names(ranks_db_hc) <- gseaDat_db_hc$Entrez
head(ranks_db_hc)
tail(ranks_db_hc)

# Plot ranks
barplot(sort(ranks_db_hc, decreasing = T))

# Load pathways
pathwaysH <- Hs.c7

# Conduct analysis
fgseaRes_db_hc <- fgsea(pathwaysH, ranks_db_hc, minSize = 15, maxSize = 500)

# Examining the top 10 results
head(fgseaRes_db_hc)

# GSEA table plot
topUp <- fgseaRes_db_hc %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes_db_hc %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway],
              ranks_db_hc, 
              fgseaRes_db_hc, 
              gseaParam = 0.1)

#supportedOrganisms() %>% filter(str_detect(Genome, "hs"))

# List of differentially expressed genes
isSigGene <- gseaDat_db_hc$padj < 0.01 & !is.na(gseaDat_db_hc$padj)
genes <- as.integer(isSigGene)
names(genes) <- gseaDat_db_hc$GeneID

# Probability weighting function
pwf <- nullp(genes, "hg19", "ensGene", bias.data = gseaDat_db_hc$medianTxLength)

# GO enrichment analysis
goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:BP"))

# Plot the top 10
g.go.db_hc <- goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "DB versus pHC") + theme_bw()

ggsave("g.go.db_hc.pdf", g.go.db_hc, dpi = 300, height = 6, width = 8)

#### DV versus HC ####
# dds_samples_new_dv_hc <- dds_samples
# dds_samples_new_dv_hc <- DESeq(dds_samples_new_dv_hc)
res_dv_hc <- results(dds_samples_new, contrast = c("phenotype", "HC", "DV"))
res_dv_hc

res_dv_hc$GeneID <- row.names(res_dv_hc)

# Create ranks
gseaDat_dv_hc <- res_dv_hc
gseaDat_dv_hc$Entrez = mapIds(org.Hs.eg.db,
                              keys=row.names(gseaDat_dv_hc),
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first") #Ensembl gene IDs to Entrez gene IDs

gseaDat_dv_hc <- gseaDat_dv_hc[!is.na(gseaDat_dv_hc$Entrez),]

gseaDat_dv_hc$Symbol = mapIds(org.Hs.eg.db,
                              keys=gseaDat_dv_hc$Entrez,
                              column="SYMBOL",
                              keytype="ENTREZID",
                              multiVals="first") #Entrez gene IDs to gene symbols

ranks_dv_hc <- gseaDat_dv_hc$log2FoldChange
names(ranks_dv_hc) <- gseaDat_dv_hc$Entrez
head(ranks_dv_hc)
tail(ranks_dv_hc)

# Plot ranks
barplot(sort(ranks_dv_hc, decreasing = T))

# Load pathways
pathwaysH <- Hs.c7

# Conduct analysis
fgseaRes_dv_hc <- fgsea(pathwaysH, ranks_dv_hc, minSize = 15, maxSize = 500)

# Examining the top 10 results
head(fgseaRes_dv_hc)

# GSEA table plot
topUp <- fgseaRes_dv_hc %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes_dv_hc %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway],
              ranks_dv_hc, 
              fgseaRes_dv_hc, 
              gseaParam = 0.1)

#supportedOrganisms() %>% filter(str_detect(Genome, "hs"))

# List of differentially expressed genes
isSigGene <- gseaDat_dv_hc$padj < 0.01 & !is.na(gseaDat_dv_hc$padj)
genes <- as.integer(isSigGene)
names(genes) <- gseaDat_dv_hc$GeneID

# Probability weighting function
pwf <- nullp(genes, "hg19", "ensGene", bias.data = gseaDat_dv_hc$medianTxLength)

# GO enrichment analysis
goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:BP"))

# Plot the top 10
g.go.dv_hc <- goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "DV versus pHC") + theme_bw()

ggsave("g.go.dv_hc.pdf", g.go.dv_hc, dpi = 300, height = 6, width = 8)

#### MIS-C versus DB ####
# dds_samples_new_dv_hc <- dds_samples
# dds_samples_new_dv_hc <- DESeq(dds_samples_new_dv_hc)
res_db_misc <- results(dds_samples_new, contrast = c("phenotype", "DB", "MIS_C"))
res_db_misc

res_db_misc$GeneID <- row.names(res_db_misc)

# Create ranks
gseaDat_db_misc <- res_db_misc
gseaDat_db_misc$Entrez = mapIds(org.Hs.eg.db,
                                keys=row.names(gseaDat_db_misc),
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first") #Ensembl gene IDs to Entrez gene IDs

gseaDat_db_misc <- gseaDat_db_misc[!is.na(gseaDat_db_misc$Entrez),]

gseaDat_db_misc$Symbol = mapIds(org.Hs.eg.db,
                                keys=gseaDat_db_misc$Entrez,
                                column="SYMBOL",
                                keytype="ENTREZID",
                                multiVals="first") #Entrez gene IDs to gene symbols

ranks_db_misc <- gseaDat_db_misc$log2FoldChange
names(ranks_db_misc) <- gseaDat_db_misc$Entrez
head(ranks_db_misc)
tail(ranks_db_misc)

# Plot ranks
barplot(sort(ranks_db_misc, decreasing = T))

# Load pathways
pathwaysH <- Hs.c7

# Conduct analysis
fgseaRes_db_misc <- fgsea(pathwaysH, ranks_db_misc, minSize = 15, maxSize = 500)

# Examining the top 10 results
head(fgseaRes_db_misc)

# GSEA table plot
topUp <- fgseaRes_db_misc %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes_db_misc %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway],
              ranks_db_misc, 
              fgseaRes_db_misc, 
              gseaParam = 0.1)

#supportedOrganisms() %>% filter(str_detect(Genome, "hs"))

# List of differentially expressed genes
isSigGene <- gseaDat_db_misc$padj < 0.01 & !is.na(gseaDat_db_misc$padj)
genes <- as.integer(isSigGene)
names(genes) <- gseaDat_db_misc$GeneID

# Probability weighting function
pwf <- nullp(genes, "hg19", "ensGene", bias.data = gseaDat_db_misc$medianTxLength)

# GO enrichment analysis
goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:BP"))

# Plot the top 10
g.go.db_misc <- goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title = "MIS-C versus DB") + theme_bw()

ggsave("g.go.db_misc.pdf", g.go.db_misc, dpi = 300, height = 6, width = 8)
