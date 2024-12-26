# Load required libraries
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorlevelPlot)
library(gridExtra)
library(data.table)
library(ComplexHeatmap)
library(circlize)

# Allow multithreading for WGCNA
allowWGCNAThreads()

# Fetch data
geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))

# Select relevant phenotype columns
phenoData <- phenoData[, c(1, 2, 46:50)]
head(phenoData)

# Load count data
data <- read.delim('/Users/sruthi/Desktop/GSE152418_p20047_Study1_RawCounts (1).txt', header = TRUE)

# Prepare data for analysis
data <- data %>%
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>%
  mutate(samples = gsub('\\.', '-', samples)) %>%
  inner_join(phenoData, by = c('samples' = 'title')) %>%
  select(1, 3, 4) %>%
  spread(key = 'geo_accession', value = 'counts') %>%
  column_to_rownames(var = 'ENSEMBLID')

# Outlier detection
gsg <- goodSamplesGenes(t(data))
if (!gsg$allOK) {
  data <- data[gsg$goodGenes, ]
}

# PCA for outlier visualization
pca <- prcomp(t(data))
pca.dat <- as.data.frame(pca$x)
pca.var <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 2)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0("PC1: ", pca.var[1], "%"), y = paste0("PC2: ", pca.var[2], "%"))

# Remove specified outliers
samples_to_be_excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
data <- data[, !colnames(data) %in% samples_to_be_excluded]

# Prepare metadata
colData <- phenoData[!rownames(phenoData) %in% samples_to_be_excluded, ]
colData <- colData %>%
  rename_with(~ gsub(":ch1", "", .), everything()) %>%
  rename_with(~ gsub("\\s", "_", .), everything())

# Ensure row and column alignment
stopifnot(all(rownames(colData) %in% colnames(data)))
stopifnot(all(rownames(colData) == colnames(data)))

# Create DESeq2 dataset
colData$condition <- factor(ifelse(colData$severity == "ICU", "Severe",
                                   ifelse(colData$severity %in% c("Moderate", "Severe"), "Moderate",
                                          "Convalescent")))
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)

# Differential expression analysis
dds <- DESeq(dds)
deResults <- results(dds, contrast = c("condition", "Severe", "Moderate"))

# Filter significant DEGs
sigResults <- deResults[which(deResults$pvalue < 0.05 & abs(deResults$log2FoldChange) > log2(2)), ]
numUpregulated <- sum(sigResults$log2FoldChange > log2(2))
numDownregulated <- sum(sigResults$log2FoldChange < -log2(2))

# Visualize DEGs
ggplot(data.frame(Category = c("Upregulated", "Downregulated"), Count = c(numUpregulated, numDownregulated)),
       aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Number of Differentially Expressed Genes", y = "Count") +
  theme_minimal()

# Normalize data for WGCNA
dds_norm <- vst(dds)
norm.counts <- assay(dds_norm)

# Choose soft threshold for WGCNA
sft <- pickSoftThreshold(norm.counts, powerVector = c(1:10, seq(12, 50, by = 2)), networkType = "signed")
ggplot(sft$fitIndices, aes(power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Soft Threshold', y = 'Scale-Free Topology Model Fit') +
  theme_classic()

# Construct WGCNA network
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = 18,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# Plot module dendrogram and colors
plotDendroAndColors(bwnet$dendrograms[[1]],
                    bwnet$colors,
                    dendroLabels = FALSE,
                    addGuide = TRUE)

# Module-trait relationships
traits <- colData %>%
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>%
  select(disease_state_bin)

module.traits.corr <- cor(bwnet$MEs, traits, use = "pairwise.complete.obs")
Heatmap(module.traits.corr,
        name = "Correlation",
        row_title = "Modules",
        column_title = "Traits",
        col = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (!is.na(module.traits.corr[i, j])) {
            grid.text(sprintf("%.2f", module.traits.corr[i, j]), x, y, gp = gpar(fontsize = 10))
          }
        })
