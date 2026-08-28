#!/usr/bin/env Rscript
library(argparser)
library(parallel)
library(BiocGenerics)
library(Biobase)
library(genefilter)
library(WGCNA)
library(flashClust)
library(caret)

PROJECT_DIR <- getwd()
EXPRESSION_FILE <- file.path(PROJECT_DIR, "merge.TPM.txt")
OUTPUT_DIR <- file.path(PROJECT_DIR, "WGCNA_output")

# Set to a trait file path to perform module-trait association analysis.
# Keep NULL if only module identification is needed.
TRAIT_FILE <- NULL
# Example:
# TRAIT_FILE <- file.path(PROJECT_DIR, "trait.txt")

# Expression filtering. IQR_QUANTILE = 0.5 retains genes whose IQR is above
# the median IQR, matching the intent of the original varFilter setting.
FILTER_BY_IQR <- TRUE
IQR_QUANTILE <- 0.5

# Sample outlier removal. NA means that samples are not removed automatically.
# After inspecting sampleClustering.pdf, replace NA with a justified cut height
# if automatic removal is required, for example: OUTLIER_CUT_HEIGHT <- 15000.
OUTLIER_CUT_HEIGHT <- NA_real_
MIN_SAMPLE_CLUSTER_SIZE <- 10

# Soft-threshold selection. NA selects the first power reaching the target R2;
# if no power reaches the target, a sample-size-based fallback is used.
SOFT_POWER <- NA_integer_
SCALE_FREE_R2_THRESHOLD <- 0.80
POWER_CANDIDATES <- 1:20

# Network construction.
NETWORK_TYPE <- "unsigned"
TOM_TYPE <- "unsigned"
MAX_BLOCK_SIZE <- 30000
MIN_MODULE_SIZE <- 50
MERGE_CUT_HEIGHT <- 0.25
DEEP_SPLIT <- 2
N_THREADS <- 12
SAVE_TOMS <- TRUE

# Module-expression heatmaps may be slow and very large when many modules or
# genes are present, so they are disabled by default.
PLOT_MODULE_EXPRESSION_HEATMAPS <- FALSE

# ---------------------------- End configuration -------------------------------

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("R package 'WGCNA' is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(WGCNA))

options(stringsAsFactors = FALSE)

if (!file.exists(EXPRESSION_FILE)) {
  stop("Expression file not found: ", EXPRESSION_FILE, call. = FALSE)
}
if (!is.null(TRAIT_FILE) && !file.exists(TRAIT_FILE)) {
  stop("Trait file not found: ", TRAIT_FILE, call. = FALSE)
}
if (!is.numeric(IQR_QUANTILE) || length(IQR_QUANTILE) != 1 ||
    !is.finite(IQR_QUANTILE) || IQR_QUANTILE < 0 || IQR_QUANTILE > 1) {
  stop("IQR_QUANTILE must be between 0 and 1.", call. = FALSE)
}
if (!NETWORK_TYPE %in% c("unsigned", "signed", "signed hybrid")) {
  stop("Unsupported NETWORK_TYPE: ", NETWORK_TYPE, call. = FALSE)
}
if (!TOM_TYPE %in% c("unsigned", "signed", "signed Nowick")) {
  stop("Unsupported TOM_TYPE: ", TOM_TYPE, call. = FALSE)
}
if (!is.numeric(MERGE_CUT_HEIGHT) || MERGE_CUT_HEIGHT <= 0 ||
    MERGE_CUT_HEIGHT >= 1) {
  stop("MERGE_CUT_HEIGHT must be between 0 and 1.", call. = FALSE)
}
if (!is.numeric(N_THREADS) || N_THREADS < 1) {
  stop("N_THREADS must be a positive integer.", call. = FALSE)
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUTPUT_DIR <- normalizePath(OUTPUT_DIR, mustWork = TRUE)
EXPRESSION_FILE <- normalizePath(EXPRESSION_FILE, mustWork = TRUE)
if (!is.null(TRAIT_FILE)) TRAIT_FILE <- normalizePath(TRAIT_FILE, mustWork = TRUE)

tryCatch(
  allowWGCNAThreads(nThreads = as.integer(N_THREADS)),
  error = function(error) {
    warning("Multithreading could not be enabled; continuing with available threads: ",
            conditionMessage(error))
  }
)

message("Reading expression matrix: ", EXPRESSION_FILE)
expression_raw <- read.delim(
  EXPRESSION_FILE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(expression_raw) < 2 || ncol(expression_raw) < 4) {
  stop("The expression matrix must contain at least 2 genes and 4 samples.",
       call. = FALSE)
}

raw_matrix <- as.matrix(expression_raw)
expression_matrix <- suppressWarnings(matrix(
  as.numeric(raw_matrix),
  nrow = nrow(raw_matrix),
  ncol = ncol(raw_matrix),
  dimnames = dimnames(raw_matrix)
))

invalid_numeric <- is.na(expression_matrix) & !is.na(raw_matrix) &
  nzchar(trimws(raw_matrix)) & toupper(trimws(raw_matrix)) != "NA"
if (any(invalid_numeric)) {
  first_bad <- which(invalid_numeric, arr.ind = TRUE)[1, ]
  stop(
    "A non-numeric expression value was found at gene '",
    rownames(raw_matrix)[first_bad[[1]]], "', sample '",
    colnames(raw_matrix)[first_bad[[2]]], "'.",
    call. = FALSE
  )
}
expression_matrix[!is.finite(expression_matrix)] <- NA_real_

# Filter low-variability genes before transposing to the WGCNA format.
if (FILTER_BY_IQR) {
  gene_iqr <- apply(expression_matrix, 1, IQR, na.rm = TRUE)
  iqr_threshold <- as.numeric(quantile(
    gene_iqr[is.finite(gene_iqr)],
    probs = IQR_QUANTILE,
    na.rm = TRUE,
    names = FALSE
  ))
  keep_by_iqr <- is.finite(gene_iqr) & gene_iqr > iqr_threshold
  if (!any(keep_by_iqr)) {
    stop("No genes remained after IQR filtering.", call. = FALSE)
  }
  filtered_matrix <- expression_matrix[keep_by_iqr, , drop = FALSE]
} else {
  gene_iqr <- apply(expression_matrix, 1, IQR, na.rm = TRUE)
  iqr_threshold <- NA_real_
  keep_by_iqr <- rep(TRUE, nrow(expression_matrix))
  filtered_matrix <- expression_matrix
}

filter_summary <- data.frame(
  Metric = c(
    "Genes_before_filtering", "Genes_after_IQR_filtering",
    "Samples_before_QC", "IQR_filter_enabled", "IQR_quantile",
    "IQR_threshold"
  ),
  Value = c(
    nrow(expression_matrix), nrow(filtered_matrix), ncol(filtered_matrix),
    FILTER_BY_IQR, IQR_QUANTILE, iqr_threshold
  ),
  stringsAsFactors = FALSE
)
write.table(
  filter_summary,
  file.path(OUTPUT_DIR, "filter_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# WGCNA requires samples in rows and genes in columns.
datExpr0 <- as.data.frame(t(filtered_matrix), check.names = FALSE)

quality <- goodSamplesGenes(datExpr0, verbose = 3)
if (!quality$allOK) {
  removed_genes <- colnames(datExpr0)[!quality$goodGenes]
  removed_samples <- rownames(datExpr0)[!quality$goodSamples]

  if (length(removed_genes) > 0) {
    writeLines(removed_genes, file.path(OUTPUT_DIR, "removed_genes.txt"))
  }
  if (length(removed_samples) > 0) {
    writeLines(removed_samples, file.path(OUTPUT_DIR, "removed_samples_QC.txt"))
  }
  datExpr0 <- datExpr0[quality$goodSamples, quality$goodGenes, drop = FALSE]
}

if (nrow(datExpr0) < 4 || ncol(datExpr0) < MIN_MODULE_SIZE) {
  stop(
    "Too few samples or genes remain after quality control: ",
    nrow(datExpr0), " samples and ", ncol(datExpr0), " genes.",
    call. = FALSE
  )
}

sampleTree <- hclust(dist(datExpr0, method = "manhattan"), method = "average")
pdf(file.path(OUTPUT_DIR, "sampleClustering.pdf"), width = 12, height = 9)
par(cex = 1, mar = c(2, 4, 2, 0))
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.2,
  cex.axis = 1.2,
  cex.main = 1.5
)
if (is.finite(OUTLIER_CUT_HEIGHT)) {
  abline(h = OUTLIER_CUT_HEIGHT, col = "red", lty = 2)
}
dev.off()

keep_samples <- rep(TRUE, nrow(datExpr0))
names(keep_samples) <- rownames(datExpr0)

if (is.finite(OUTLIER_CUT_HEIGHT)) {
  clusters <- cutreeStatic(
    sampleTree,
    cutHeight = OUTLIER_CUT_HEIGHT,
    minSize = MIN_SAMPLE_CLUSTER_SIZE
  )
  valid_clusters <- clusters[clusters > 0]
  if (length(valid_clusters) == 0) {
    warning("No valid sample cluster was found; all QC-passing samples were retained.")
  } else {
    largest_cluster <- as.integer(names(which.max(table(valid_clusters))))
    keep_samples <- clusters == largest_cluster
    removed_outliers <- rownames(datExpr0)[!keep_samples]
    if (length(removed_outliers) > 0) {
      writeLines(
        removed_outliers,
        file.path(OUTPUT_DIR, "removed_samples_clustering.txt")
      )
    }
  }
}

datExpr <- datExpr0[keep_samples, , drop = FALSE]
nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

if (nSamples < 4 || nGenes < MIN_MODULE_SIZE) {
  stop(
    "Too few samples or genes remain after outlier removal: ",
    nSamples, " samples and ", nGenes, " genes.",
    call. = FALSE
  )
}

write.table(
  t(datExpr),
  file.path(OUTPUT_DIR, "InputData_processed.tsv"),
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

message("Selecting soft-threshold power...")
sft <- pickSoftThreshold(
  datExpr,
  powerVector = POWER_CANDIDATES,
  networkType = NETWORK_TYPE,
  verbose = 5
)

fit_indices <- sft$fitIndices
signed_r2 <- -sign(fit_indices[, 3]) * fit_indices[, 2]

recommended_power <- function(sample_number) {
  if (sample_number < 20) return(10L)
  if (sample_number < 30) return(9L)
  if (sample_number < 40) return(8L)
  if (sample_number < 60) return(7L)
  6L
}

if (length(SOFT_POWER) == 1 && is.finite(SOFT_POWER)) {
  powerEstimate <- as.integer(SOFT_POWER)
  if (powerEstimate < 1) stop("SOFT_POWER must be a positive integer or NA.")
} else {
  passing_indices <- which(
    is.finite(signed_r2) & signed_r2 >= SCALE_FREE_R2_THRESHOLD
  )
  if (length(passing_indices) > 0) {
    powerEstimate <- as.integer(fit_indices[passing_indices[[1]], 1])
  } else {
    powerEstimate <- recommended_power(nSamples)
    warning(
      "No tested power reached the target scale-free R2; using power ",
      powerEstimate, " based on sample size."
    )
  }
}

soft_threshold_table <- data.frame(
  fit_indices,
  Signed_R2 = signed_r2,
  Selected = fit_indices[, 1] == powerEstimate,
  check.names = FALSE
)
write.table(
  soft_threshold_table,
  file.path(OUTPUT_DIR, "soft_threshold_statistics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(file.path(OUTPUT_DIR, "soft-thresholding.power.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(
  fit_indices[, 1], signed_r2,
  xlab = "Soft-threshold power",
  ylab = "Scale-free topology model fit, signed R2",
  type = "n",
  main = "Scale independence"
)
text(fit_indices[, 1], signed_r2, labels = fit_indices[, 1], cex = 0.9, col = "red")
abline(h = SCALE_FREE_R2_THRESHOLD, col = "red", lty = 2)
abline(v = powerEstimate, col = "blue", lty = 2)
plot(
  fit_indices[, 1], fit_indices[, 5],
  xlab = "Soft-threshold power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  fit_indices[, 1], fit_indices[, 5],
  labels = fit_indices[, 1], cex = 0.9, col = "red"
)
abline(v = powerEstimate, col = "blue", lty = 2)
dev.off()

message("Constructing network with power = ", powerEstimate, "...")
tom_file_base <- file.path(OUTPUT_DIR, "blockTOM")
net <- blockwiseModules(
  datExpr,
  maxBlockSize = MAX_BLOCK_SIZE,
  power = powerEstimate,
  networkType = NETWORK_TYPE,
  TOMType = TOM_TYPE,
  saveTOMs = SAVE_TOMS,
  saveTOMFileBase = tom_file_base,
  minModuleSize = MIN_MODULE_SIZE,
  deepSplit = DEEP_SPLIT,
  reassignThreshold = 0,
  mergeCutHeight = MERGE_CUT_HEIGHT,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  nThreads = as.integer(N_THREADS),
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# Plot every block, rather than only the first block.
pdf(file.path(OUTPUT_DIR, "Gene_ClusterDendrogram.pdf"), width = 12, height = 9)
for (block_index in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[block_index]]
  plotDendroAndColors(
    net$dendrograms[[block_index]],
    moduleColors[block_genes],
    paste0("Module colors - block ", block_index),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
}
dev.off()

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

geneModuleMembership <- as.data.frame(
  cor(datExpr, MEs, use = "pairwise.complete.obs")
)
MMPvalue <- as.data.frame(
  corPvalueStudent(as.matrix(geneModuleMembership), nSamples)
)
names(geneModuleMembership) <- sub("^ME", "MM_", names(geneModuleMembership))
names(MMPvalue) <- sub("^ME", "P_MM_", names(MMPvalue))

geneModule <- data.frame(
  Gene_ID = colnames(datExpr),
  Module_Label = moduleLabels,
  Module_Color = moduleColors,
  geneModuleMembership,
  MMPvalue,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(
  geneModule,
  file.path(OUTPUT_DIR, "gene_module_membership.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

module_sizes <- as.data.frame(table(Module_Color = moduleColors))
names(module_sizes)[[2]] <- "Gene_Count"
write.table(
  module_sizes,
  file.path(OUTPUT_DIR, "module_sizes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  MEs,
  file.path(OUTPUT_DIR, "module_eigengenes.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

if (ncol(MEs) >= 2) {
  pdf(file.path(OUTPUT_DIR, "Eigengene_Dendrogram.pdf"), width = 7, height = 7)
  plotEigengeneNetworks(
    MEs,
    "",
    marDendro = c(0, 4, 2, 2),
    marHeatmap = c(3, 4, 2, 2),
    xLabelsAngle = 90
  )
  dev.off()
}

# Optional module-trait association analysis.
moduleTraitCor <- NULL
moduleTraitPvalue <- NULL
if (!is.null(TRAIT_FILE)) {
  message("Reading trait matrix: ", TRAIT_FILE)
  trait_raw <- read.delim(
    TRAIT_FILE,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  trait_matrix_raw <- as.matrix(trait_raw)
  trait_matrix <- suppressWarnings(matrix(
    as.numeric(trait_matrix_raw),
    nrow = nrow(trait_matrix_raw),
    ncol = ncol(trait_matrix_raw),
    dimnames = dimnames(trait_matrix_raw)
  ))
  invalid_traits <- is.na(trait_matrix) & !is.na(trait_matrix_raw) &
    nzchar(trimws(trait_matrix_raw)) & toupper(trimws(trait_matrix_raw)) != "NA"
  if (any(invalid_traits)) {
    stop("The trait file contains non-numeric values.", call. = FALSE)
  }

  trait_data <- as.data.frame(trait_matrix, check.names = FALSE)
  common_samples <- intersect(rownames(datExpr), rownames(trait_data))
  if (length(common_samples) < 4) {
    stop("Fewer than four samples are shared by expression and trait data.",
         call. = FALSE)
  }

  trait_data <- trait_data[common_samples, , drop = FALSE]
  trait_variance <- vapply(trait_data, var, numeric(1), na.rm = TRUE)
  valid_traits <- is.finite(trait_variance) & trait_variance > 0
  if (!any(valid_traits)) {
    stop("No variable numeric traits remain for module-trait analysis.",
         call. = FALSE)
  }
  trait_data <- trait_data[, valid_traits, drop = FALSE]
  me_for_traits <- MEs[common_samples, , drop = FALSE]

  complete_rows <- complete.cases(trait_data) & complete.cases(me_for_traits)
  trait_data <- trait_data[complete_rows, , drop = FALSE]
  me_for_traits <- me_for_traits[complete_rows, , drop = FALSE]
  trait_sample_number <- nrow(trait_data)
  if (trait_sample_number < 4) {
    stop("Fewer than four complete samples remain for module-trait analysis.",
         call. = FALSE)
  }

  moduleTraitCor <- cor(me_for_traits, trait_data, use = "everything")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, trait_sample_number)

  association_grid <- expand.grid(
    Module = rownames(moduleTraitCor),
    Trait = colnames(moduleTraitCor),
    stringsAsFactors = FALSE
  )
  association_grid$Correlation <- as.vector(moduleTraitCor)
  association_grid$P_value <- as.vector(moduleTraitPvalue)
  association_grid$FDR_BH <- p.adjust(association_grid$P_value, method = "BH")
  association_grid$N <- trait_sample_number
  write.table(
    association_grid,
    file.path(OUTPUT_DIR, "module_trait_associations.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  text_matrix <- paste0(
    signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")"
  )
  dim(text_matrix) <- dim(moduleTraitCor)

  heatmap_width <- max(7, min(18, 4 + 0.65 * ncol(moduleTraitCor)))
  heatmap_height <- max(6, min(18, 3 + 0.35 * nrow(moduleTraitCor)))
  pdf(
    file.path(OUTPUT_DIR, "module_trait_relationships.pdf"),
    width = heatmap_width,
    height = heatmap_height
  )
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(moduleTraitCor),
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    xLabelsAngle = 45,
    colorLabels = FALSE,
    colors = blueWhiteRed(100),
    textMatrix = text_matrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    zlim = c(-1, 1),
    main = "Module-trait relationships"
  )
  dev.off()
}

if (PLOT_MODULE_EXPRESSION_HEATMAPS) {
  modules_to_plot <- setdiff(unique(moduleColors), "grey")
  for (module in modules_to_plot) {
    module_genes <- moduleColors == module
    module_me_name <- paste0("ME", module)
    if (!module_me_name %in% names(MEs)) next

    output_file <- file.path(
      OUTPUT_DIR,
      paste0("Heatmap_EigengeneExpression_", module, ".pdf")
    )
    pdf(output_file, width = 10, height = 8)
    par(mfrow = c(2, 1), mar = c(0.3, 5.5, 3, 2))
    plotMat(
      t(scale(datExpr[, module_genes, drop = FALSE])),
      nrgcols = 30,
      rlabels = FALSE,
      clabels = FALSE,
      rcols = module,
      main = module,
      cex.main = 1.5
    )
    par(mar = c(5, 4.2, 0, 0.7))
    module_me <- MEs[, module_me_name]
    names(module_me) <- rownames(datExpr)
    barplot(
      module_me,
      col = module,
      ylab = "Eigengene expression",
      xlab = "Sample"
    )
    dev.off()
  }
}

run_summary <- data.frame(
  Metric = c(
    "Expression_file", "Trait_file", "Samples_final", "Genes_final",
    "Selected_power", "Network_type", "TOM_type", "Module_count_excluding_grey"
  ),
  Value = c(
    EXPRESSION_FILE,
    if (is.null(TRAIT_FILE)) "Not supplied" else TRAIT_FILE,
    nSamples,
    nGenes,
    powerEstimate,
    NETWORK_TYPE,
    TOM_TYPE,
    length(setdiff(unique(moduleColors), "grey"))
  ),
  stringsAsFactors = FALSE
)
write.table(
  run_summary,
  file.path(OUTPUT_DIR, "run_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

save(
  datExpr,
  net,
  MEs,
  geneModule,
  moduleLabels,
  moduleColors,
  geneModuleMembership,
  MMPvalue,
  moduleTraitCor,
  moduleTraitPvalue,
  file = file.path(OUTPUT_DIR, "WGCNA_Module_Result.RData")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)

message("WGCNA analysis completed successfully.")
message("Results directory: ", OUTPUT_DIR)
