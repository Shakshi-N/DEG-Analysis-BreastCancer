## =============================================================================
## Project Title: Differential Gene Expression Analysis of Breast Tumor vs Normal
##                Tissue using RNA-seq (GSE183947)
##
## Pipeline: Data acquisition -> Metadata reshaping -> Count matrix construction
##           -> DESeq2 DGE analysis -> Visualization -> Functional enrichment
##           -> ML-based classification (Random Forest) on top DEGs
## =============================================================================

## -----------------------------------------------------------------------------
## 0. PACKAGES
## -----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "org.Hs.eg.db", "clusterProfiler",
                        "enrichplot", "GEOquery", "EnhancedVolcano"),
                      update = FALSE, ask = FALSE)

install.packages(c("tidyverse", "pheatmap", "enrichR", "randomForest",
                    "caret", "e1071"), dependencies = TRUE)

library(tidyverse)      # data wrangling + ggplot2
library(GEOquery)       # metadata retrieval from GEO
library(DESeq2)         # differential expression
library(EnhancedVolcano)# volcano plots
library(pheatmap)       # heatmaps
library(org.Hs.eg.db)   # gene ID mapping
library(clusterProfiler)# GO / functional enrichment
library(enrichR)        # GO / KEGG enrichment (web API)
library(randomForest)   # ML classifier
library(caret)          # train/test split + confusion matrix

set.seed(123)

## Ensure output directories exist regardless of how the repo was cloned
## (git does not track empty folders, so these are recreated on first run)
dir.create("results/plots",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data",           recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. DATA ACQUISITION AND PREPARATION
## -----------------------------------------------------------------------------
## Dataset: GSE183947 (NCBI GEO) - 60 breast tissue samples
##          30 breast tumor + 30 matched normal tissue
##          Already pre-normalized in FPKM units by original submitters

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

# FPKM expression matrix (genes x samples), downloaded from GEO supplementary files
dat <- read.csv(file = "data/GSE183947_fpkm.csv")

# Series matrix / phenotype metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

## -----------------------------------------------------------------------------
## 2. METADATA PROCESSING AND RESHAPING
## -----------------------------------------------------------------------------
# Keep only relevant columns and rename for clarity (tissue & metastasis status)
metadata.modified <- metadata %>%
  select(1, 10, 11, 17) %>%
  rename(tissue     = characteristics_ch1,
         metastasis = characteristics_ch1.1) %>%
  mutate(tissue     = gsub("tissue: ", "", tissue),
         metastasis = gsub("metastasis: ", "", metastasis))

# Reshape wide FPKM matrix -> long format (one row per gene-sample pair)
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "samples", value = "FPKM", -gene)

# Merge expression data with cleaned metadata using sample identifiers
dat.long <- dat.long %>%
  left_join(metadata.modified, by = c("samples" = "description"))

## -----------------------------------------------------------------------------
## 3. COUNT MATRIX & SAMPLE INFORMATION
## -----------------------------------------------------------------------------
# DESeq2 requires integer counts. FPKM values are rounded to nearest integer
# to approximate a count matrix (as done in the original assignment).
count_matrix <- dat.long %>%
  select(gene, samples, FPKM) %>%
  pivot_wider(names_from = samples, values_from = FPKM) %>%
  remove_rownames() %>%
  column_to_rownames("gene") %>%
  as.matrix()

count_matrix_int <- round(count_matrix)

# Sample info table: assign each sample to condition (Tumor / Normal)
sample_info <- metadata.modified %>%
  select(description, tissue, metastasis) %>%
  distinct() %>%
  remove_rownames() %>%                # metadata.modified inherits GSM
  column_to_rownames("description") %>% # accession rownames from pData();
  mutate(condition = ifelse(tissue == "breast tumor", "Tumor",  # strip them
                      ifelse(tissue == "normal breast tissue", "Normal", NA)))# before setting new ones

sample_info$condition <- factor(sample_info$condition, levels = c("Normal", "Tumor"))

# IMPORTANT: ensure column order of count matrix matches row order of sample_info
# (DESeq2 requires this alignment; do NOT assume the original ordering holds)
stopifnot(all(colnames(count_matrix_int) %in% rownames(sample_info)))
sample_info <- sample_info[colnames(count_matrix_int), , drop = FALSE]
stopifnot(all(colnames(count_matrix_int) == rownames(sample_info)))

## -----------------------------------------------------------------------------
## 4. DIFFERENTIAL GENE EXPRESSION ANALYSIS (DESeq2)
## -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
  colData   = sample_info,
  design    = ~ condition
)

# Filter low-count genes (reduce noise / multiple-testing burden)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Normalize + fit negative binomial GLM per gene
dds <- DESeq(dds)

# Extract results: Tumor vs Normal
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- res[order(res$padj), ]

# Significant DEGs: padj < 0.05 & |log2FC| >= 1
res_sig   <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)
up_genes   <- rownames(subset(res_sig, log2FoldChange > 1))
down_genes <- rownames(subset(res_sig, log2FoldChange < -1))

write.csv(as.data.frame(res),     "results/tables/DESeq2_full_results.csv")
write.csv(as.data.frame(res_sig), "results/tables/DESeq2_significant_DEGs.csv")

## -----------------------------------------------------------------------------
## 5. VISUALIZATION
## -----------------------------------------------------------------------------

## 5a. Boxplots of top 5 DEGs (by padj)
top5_genes <- res_sig %>%
  as.data.frame() %>%
  arrange(padj) %>%
  head(5) %>%
  rownames()

norm_counts <- counts(dds, normalized = TRUE)

top5_counts <- norm_counts[top5_genes, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count")

top5_counts <- top5_counts %>%
  left_join(sample_info %>% rownames_to_column("sample"), by = "sample")

p_box <- ggplot(top5_counts, aes(x = condition, y = log2(count + 1), fill = condition)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  labs(y = "Log2 Normalized Counts", x = "Condition",
       title = "Top 5 DEGs (Tumor vs Normal)") +
  scale_fill_manual(values = c("Normal" = "skyblue", "Tumor" = "salmon"))

ggsave("results/plots/top5_DEGs_boxplot.png", p_box, width = 9, height = 6, dpi = 300)

## 5b. MA plot
png("results/plots/MA_plot.png", width = 1600, height = 1200, res = 200)
plotMA(res, ylim = c(-5, 5))
dev.off()

## 5c. Volcano plot
png("results/plots/volcano_plot.png", width = 2000, height = 1600, res = 220)
print(
  EnhancedVolcano(res,
                   lab = rownames(res),
                   x = "log2FoldChange",
                   y = "padj",
                   xlab = "Log2 Fold Change",
                   ylab = "-Log10 adjusted p-value",
                   pCutoff = 0.05,
                   FCcutoff = 2,
                   pointSize = 2.5,
                   labSize = 3.5,
                   colAlpha = 0.6,
                   legendLabels = c("NS", "Log2FC", "p-adj", "p-adj & Log2FC"),
                   legendPosition = "right",
                   title = "Tumor vs Normal (DEGs)",
                   subtitle = "DESeq2 analysis")
)
dev.off()

## 5d. Heatmap of top 50 DEGs
top50_genes <- rownames(res_sig)[1:min(50, nrow(res_sig))]
mat <- norm_counts[top50_genes, , drop = FALSE]
mat <- log2(mat + 1)
annotation <- sample_info[, "condition", drop = FALSE]

png("results/plots/heatmap_top50_DEGs.png", width = 2000, height = 2400, res = 220)
pheatmap(mat,
         annotation_col = annotation,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row")
dev.off()

## -----------------------------------------------------------------------------
## 6. FUNCTIONAL ENRICHMENT ANALYSIS
## -----------------------------------------------------------------------------

## 6a. GO Biological Process enrichment (clusterProfiler) - upregulated genes
up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")

enrich_up <- enrichGO(gene          = up_entrez$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "BP",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2,
                       pAdjustMethod = "BH",
                       readable      = TRUE)

# Always save the result table (useful even if 0 terms pass cutoffs)
write.csv(as.data.frame(enrich_up), "results/tables/GO_BP_enrichment_upregulated.csv",
          row.names = FALSE)

# NOTE: enrichplot::dotplot() has a long-standing, widely-reported version-
# compatibility bug with newer ggplot2 releases that throws exactly this
# error ("unique.default(x, nmax = nmax): unique() applies only to
# vectors") even when the enrichResult has valid, non-empty terms. Rather
# than depend on dotplot() internals, the plot is built directly and
# reliably from the enrichment result data frame.
n_go_terms <- nrow(as.data.frame(enrich_up))

if (n_go_terms > 0) {
  go_df <- as.data.frame(enrich_up) %>%
    arrange(p.adjust) %>%
    head(min(20, n_go_terms)) %>%
    mutate(GeneRatioNum = sapply(GeneRatio, function(r) eval(parse(text = r))))

  p_go <- ggplot(go_df, aes(x = GeneRatioNum,
                             y = reorder(Description, GeneRatioNum))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Gene Ratio", y = NULL, size = "Count", color = "p.adjust",
         title = "GO enrichment (Biological Process)") +
    theme_bw()

  ggsave("results/plots/GO_enrichment_dotplot.png", p_go,
         width = 9, height = 7, dpi = 300)
} else {
  message("No significant GO Biological Process terms found for upregulated genes ",
          "(", nrow(up_entrez), "/", length(up_genes), " genes mapped to ENTREZID). ",
          "Skipping GO dotplot - see results/tables/GO_BP_enrichment_upregulated.csv.")
}

## 6b. KEGG pathway enrichment (enrichR) - upregulated genes
dbs <- listEnrichrDbs()
enrich_dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")

enrich_up_r   <- enrichr(up_genes, enrich_dbs)
enrich_down_r <- enrichr(down_genes, enrich_dbs)

kegg_top10 <- enrich_up_r[["KEGG_2021_Human"]] %>% head(10)

p_kegg <- ggplot(kegg_top10, aes(x = reorder(Term, Combined.Score),
                                  y = Combined.Score, fill = Combined.Score)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 KEGG Pathways - Upregulated Genes",
       x = "Pathway", y = "Combined Score") +
  theme_minimal()

ggsave("results/plots/KEGG_top10_upregulated.png", p_kegg, width = 8, height = 5, dpi = 300)

## -----------------------------------------------------------------------------
## 7. MACHINE LEARNING - RANDOM FOREST CLASSIFICATION ON TOP DEGs (NESTED DESIGN)
## -----------------------------------------------------------------------------
## To evaluate the discriminative power of the identified DEGs, a Random
## Forest classifier is trained using a nested design that keeps DEG
## selection and test-set evaluation fully separated, so that no information
## from held-out samples influences which genes are chosen or how expression
## values are normalized:
##   a. Samples are split into train/test FIRST, before any DEG selection.
##   b. DESeq2 is run for DEG discovery using ONLY the training samples, to
##      pick the top 10 genes used as ML features.
##   c. Features are DESeq2-normalized, log2-transformed counts (not raw
##      counts), joined to samples by sample ID (rownames) rather than row
##      position, to avoid silent sample/label misalignment.
##   d. Test samples are normalized using size factors computed against the
##      training set's reference (via geometric means of training counts),
##      so no information from the test set - not even its own
##      normalization - enters feature selection or scaling.
## The test set is therefore genuinely held out from every step that touches
## gene/feature choice. Note that the exploratory DESeq2 analysis in Section 4
## above uses the full 60-sample dataset, as is standard for DEG discovery and
## functional enrichment reporting; this nested, sample-split design applies
## specifically to the machine learning evaluation below.

## 7a. Sample-level train/test split (before any DEG selection)
set.seed(123)
ml_trainIndex <- createDataPartition(sample_info$condition, p = 0.7, list = FALSE)
train_samples <- rownames(sample_info)[ml_trainIndex]
test_samples  <- rownames(sample_info)[-ml_trainIndex]

train_info <- sample_info[train_samples, , drop = FALSE]
test_info  <- sample_info[test_samples, , drop = FALSE]

train_counts_raw <- count_matrix_int[, train_samples, drop = FALSE]
test_counts_raw  <- count_matrix_int[, test_samples, drop = FALSE]

## 7b. DESeq2 DEG discovery using only the training samples
dds_train <- DESeqDataSetFromMatrix(
  countData = train_counts_raw,
  colData   = train_info,
  design    = ~ condition
)
dds_train <- dds_train[rowSums(counts(dds_train)) >= 10, ]
dds_train <- DESeq(dds_train)

res_train <- results(dds_train, contrast = c("condition", "Tumor", "Normal"))
res_train <- res_train[order(res_train$padj), ]
res_train_sig <- subset(res_train, padj < 0.05 & abs(log2FoldChange) >= 1)

n_available <- min(10, nrow(res_train_sig))
top10_genes_ml <- res_train_sig %>%
  as.data.frame() %>%
  arrange(padj) %>%
  slice(1:n_available) %>%
  rownames()

write.csv(as.data.frame(res_train_sig),
          "results/tables/DESeq2_train_only_DEGs_for_ML.csv")

## 7c. Training-set features: normalized, log2-transformed counts (train DESeq2 fit)
norm_train_counts <- counts(dds_train, normalized = TRUE)
train_expr <- t(log2(norm_train_counts[top10_genes_ml, , drop = FALSE] + 1)) %>%
  as.data.frame()
train_expr$condition <- train_info[rownames(train_expr), "condition"]

## 7d. Test-set features: normalize test samples AGAINST the training reference
## (geometric means from filtered training counts), so test data never
## influences its own normalization or the choice of genes.
filtered_genes   <- rownames(dds_train)
train_raw_filtered <- counts(dds_train)  # raw counts, filtered gene set, train samples

geoMeans_train <- apply(train_raw_filtered, 1, function(row) {
  if (all(row == 0)) return(0)
  exp(mean(log(row[row > 0])))
})

test_raw_filtered  <- test_counts_raw[filtered_genes, , drop = FALSE]
sizeFactors_test   <- estimateSizeFactorsForMatrix(test_raw_filtered, geoMeans = geoMeans_train)
norm_test_counts   <- sweep(test_raw_filtered, 2, sizeFactors_test, "/")

test_expr <- t(log2(norm_test_counts[top10_genes_ml, , drop = FALSE] + 1)) %>%
  as.data.frame()
test_expr$condition <- test_info[rownames(test_expr), "condition"]

## 7e. Train Random Forest classifier on training samples/features only
set.seed(123)
rf_model <- randomForest(condition ~ ., data = train_expr,
                          importance = TRUE, ntree = 500)
print(rf_model)

## 7f. Evaluate on genuinely held-out test set
pred <- predict(rf_model, newdata = test_expr)
cm <- confusionMatrix(pred, test_expr$condition)
print(cm)

## 7g. Feature importance
imp_df <- data.frame(
  Gene       = rownames(importance(rf_model)),
  Importance = importance(rf_model)[, "MeanDecreaseGini"]
)

p_imp <- ggplot(imp_df, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Random Forest Feature Importance (Train-Selected DEGs)",
       x = "Gene", y = "Mean Decrease in Gini") +
  theme_minimal()

ggsave("results/plots/RF_feature_importance.png", p_imp, width = 7, height = 5, dpi = 300)

## Save model outputs
saveRDS(rf_model, "results/tables/rf_model.rds")
write.csv(imp_df, "results/tables/RF_feature_importance.csv", row.names = FALSE)
sink("results/tables/RF_confusion_matrix.txt")
print(cm)
sink()

## -----------------------------------------------------------------------------
## 7h. ROBUSTNESS CHECK: REPEATED K-FOLD CROSS-VALIDATION (NESTED, LEAK-FREE)
## -----------------------------------------------------------------------------
## The single 70/30 split above gives one honest, leakage-free accuracy
## estimate, but with only ~18 test samples a single split carries high
## variance - a single misclassification would already move accuracy from
## 100% to ~94%. To get a more robust picture, the identical nested
## procedure (train/test split -> DEG selection on train only -> train-
## referenced test normalization -> RF) is repeated across many stratified
## folds, re-selecting DEGs from scratch within each fold's training
## partition, so no fold's evaluation is contaminated by its own test data.

cv_k       <- 5   # number of folds per repeat
cv_repeats <- 5   # number of repeats -> cv_k * cv_repeats = 25 fold fits total

set.seed(123)
cv_folds <- createMultiFolds(sample_info$condition, k = cv_k, times = cv_repeats)

run_nested_fold <- function(train_idx) {
  fold_train_samples <- rownames(sample_info)[train_idx]
  fold_test_samples  <- setdiff(rownames(sample_info), fold_train_samples)

  fold_train_info <- sample_info[fold_train_samples, , drop = FALSE]
  fold_test_info  <- sample_info[fold_test_samples, , drop = FALSE]

  fold_train_raw <- count_matrix_int[, fold_train_samples, drop = FALSE]
  fold_test_raw  <- count_matrix_int[, fold_test_samples, drop = FALSE]

  fold_dds <- DESeqDataSetFromMatrix(fold_train_raw, fold_train_info, design = ~ condition)
  fold_dds <- fold_dds[rowSums(counts(fold_dds)) >= 10, ]
  fold_dds <- suppressMessages(DESeq(fold_dds, quiet = TRUE))

  fold_res <- results(fold_dds, contrast = c("condition", "Tumor", "Normal"))
  fold_sig <- subset(fold_res, padj < 0.05 & abs(log2FoldChange) >= 1)

  n_feat <- min(10, nrow(fold_sig))
  if (n_feat < 2) return(NULL)  # not enough signal in this fold's train split - skip

  fold_top_genes <- fold_sig %>%
    as.data.frame() %>%
    arrange(padj) %>%
    slice(1:n_feat) %>%
    rownames()

  fold_norm_train <- counts(fold_dds, normalized = TRUE)
  fold_train_expr <- t(log2(fold_norm_train[fold_top_genes, , drop = FALSE] + 1)) %>%
    as.data.frame()
  fold_train_expr$condition <- fold_train_info[rownames(fold_train_expr), "condition"]

  # Normalize this fold's test samples against THIS fold's training reference only
  fold_filtered_genes <- rownames(fold_dds)
  fold_train_raw_filtered <- counts(fold_dds)
  fold_geoMeans <- apply(fold_train_raw_filtered, 1, function(row) {
    if (all(row == 0)) return(0)
    exp(mean(log(row[row > 0])))
  })

  fold_test_filtered <- fold_test_raw[fold_filtered_genes, , drop = FALSE]
  fold_test_sf     <- estimateSizeFactorsForMatrix(fold_test_filtered, geoMeans = fold_geoMeans)
  fold_norm_test   <- sweep(fold_test_filtered, 2, fold_test_sf, "/")

  fold_test_expr <- t(log2(fold_norm_test[fold_top_genes, , drop = FALSE] + 1)) %>%
    as.data.frame()
  fold_test_expr$condition <- fold_test_info[rownames(fold_test_expr), "condition"]

  fold_rf   <- randomForest(condition ~ ., data = fold_train_expr, ntree = 500)
  fold_pred <- predict(fold_rf, newdata = fold_test_expr)
  fold_cm   <- confusionMatrix(fold_pred, fold_test_expr$condition)

  data.frame(
    Accuracy    = unname(fold_cm$overall["Accuracy"]),
    Sensitivity = unname(fold_cm$byClass["Sensitivity"]),
    Specificity = unname(fold_cm$byClass["Specificity"]),
    n_features  = n_feat
  )
}

cv_results_list <- lapply(cv_folds, run_nested_fold)
cv_results <- bind_rows(cv_results_list, .id = "fold") %>% filter(!is.na(Accuracy))

write.csv(cv_results, "results/tables/RF_nested_CV_results.csv", row.names = FALSE)

cv_summary <- cv_results %>%
  summarise(
    n_folds          = n(),
    mean_accuracy    = mean(Accuracy),
    sd_accuracy      = sd(Accuracy),
    mean_sensitivity = mean(Sensitivity),
    mean_specificity = mean(Specificity)
  )
print(cv_summary)
write.csv(cv_summary, "results/tables/RF_nested_CV_summary.csv", row.names = FALSE)

p_cv <- ggplot(cv_results, aes(x = "", y = Accuracy)) +
  geom_boxplot(fill = "steelblue", width = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.05, alpha = 0.6) +
  ylim(0, 1) +
  labs(title = paste0("Nested ", cv_k, "-fold CV (x", cv_repeats, " repeats) Test Accuracy"),
       x = NULL, y = "Test Accuracy") +
  theme_minimal()

ggsave("results/plots/RF_nested_CV_accuracy.png", p_cv, width = 5, height = 5, dpi = 300)

## -----------------------------------------------------------------------------
## END OF SCRIPT
## -----------------------------------------------------------------------------
