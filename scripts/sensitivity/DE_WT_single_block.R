# Sensitivity analysis: model dependence among the parental KOLF2.2J WT replicates.
#
# The primary DE notebooks (scripts/02_DE_testing_*.ipynb) fit the JAX and MSK
# contrasts with limma-voom + duplicateCorrelation, blocking KO samples by clone
# but giving each WT replicate its own block (WT_1, WT_2, WT_3), i.e. the three
# WT replicates of the single parental line are treated as independent.
#
# Alternative (this script, default --wt-block single): WT replicates from one
# block. All WT replicates share ONE block, so the consensus within-clone
# correlation (estimated jointly from the KO clones and the WT line) also applies
# to WT. This assumes clone-to-clone variance is the same for the parental line
# and the KO clones. The consensus correlation is floored at 0: a negative
# within-clone correlation is not a valid variance component and would make
# block means look more precise than independent samples (mixed models such as
# lme4/dream constrain the clone variance to be >= 0 in the same way).
#
# --wt-block independent reproduces the original notebook model (no floor) and
# is used to verify that this script matches the exported DE tables.
#
# Everything else (sample selection, q75 gene universe, TMM, two-pass voom,
# eBayes, Storey q-values, CPM columns, DE threshold, annotation, column order)
# follows the primary notebooks.
#
# Usage:
#   Rscript DE_WT_single_block.R --out-dir <dir>                          # Alternative
#   Rscript DE_WT_single_block.R --wt-block independent --out-dir <dir>   # reproduce original
# Optional: --min-correlation <value> overrides the floor (default 0 for the
# Alternative, none for the original model).
#
# Outputs: <contrast>.tsv (same columns as the exported DE tables),
# model_summary.tsv (blocks, estimated and used correlations, DE counts),
# sessionInfo.txt.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(qvalue)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1]
}
wt_block_mode <- match.arg(get_arg("--wt-block", "single"), c("single", "independent"))
min_cor <- as.numeric(get_arg("--min-correlation",
                              if (wt_block_mode == "single") "0" else "-Inf"))
out_dir <- get_arg("--out-dir")
stopifnot(!is.null(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_dir <- "/Users/wsun/research/_by_groups/Siqi_Shen/MorPhiC"
dracc_dir <- "/Users/wsun/research/MorPhiC/exchange_experiment/DRACC-processed"
anno_file <- "/Users/wsun/research/GitHub/MorPhiC_bulk_RNAseq/data/gencode_v44_primary_assembly_info.tsv"  # same symbols as the exported tables

gene_anno <- fread(anno_file, data.table = FALSE)

# ---------------------------------------------------------------------------
# Helpers copied from the primary notebooks
# ---------------------------------------------------------------------------
fit_limma_with_optional_duplicate_correlation <- function(y, design, block = NULL, fit_label = "") {
  first_voom <- voom(y, design)
  block <- droplevels(factor(block))
  stopifnot(length(block) == ncol(first_voom$E), !anyNA(block))
  repeated_blocks <- sum(table(block) > 1)
  if (nlevels(block) > 1 && repeated_blocks > 1) {
    corfit_first <- duplicateCorrelation(first_voom, design, block = block)
    second_voom <- voom(y, design, block = block, correlation = max(corfit_first$consensus, min_cor))
    corfit_second <- duplicateCorrelation(second_voom, design, block = block)
    fit <- lmFit(second_voom, design, block = block, correlation = max(corfit_second$consensus, min_cor))
    fit$duplicate_correlation_used <- max(corfit_second$consensus, min_cor)
    fit$duplicate_correlation_first_pass <- corfit_first$consensus
    fit$duplicate_correlation <- corfit_second$consensus
    message(sprintf("  %s: duplicate-correlation (first/second pass) %.3f / %.3f",
                    fit_label, corfit_first$consensus, corfit_second$consensus))
  } else {
    fit <- lmFit(first_voom, design)
    fit$duplicate_correlation_first_pass <- NA_real_
    fit$duplicate_correlation <- NA_real_
    fit$duplicate_correlation_used <- 0
    message(sprintf("  %s: skipping duplicate-correlation (%d repeated block(s))",
                    fit_label, repeated_blocks))
  }
  fit
}

add_qvalue_and_cpm <- function(tt, counts, group, gene_info, min_cpm = 1) {
  y_cpm <- suppressMessages(calcNormFactors(DGEList(counts = counts)))
  cpm_mat <- cpm(y_cpm, log = FALSE)
  group <- droplevels(factor(group))
  group_mean_cpm <- sapply(levels(group), function(l) rowMeans(cpm_mat[, group == l, drop = FALSE]))
  colnames(group_mean_cpm) <- paste0("mean_CPM_", make.names(colnames(group_mean_cpm)))
  for (m in colnames(group_mean_cpm)) {
    tt[[m]] <- group_mean_cpm[match(rownames(tt), rownames(group_mean_cpm)), m]
  }
  cpm_filter_value <- apply(group_mean_cpm, 1, max)
  tt$max_mean_CPM <- cpm_filter_value[match(rownames(tt), names(cpm_filter_value))]
  tt$cpm_filter_value <- tt$max_mean_CPM

  p <- tt$P.Value
  qf <- qvalue(p)
  tt$qvalue <- qf$qvalues
  tt$local_fdr <- qf$lfdr
  tt$pi0_storey <- qf$pi0
  tt$pi0_lambda_0_5 <- min(1, max(0, 2 * mean(p > 0.5)))

  tt$q75 <- gene_info$q75[match(rownames(tt), rownames(gene_info))]
  tt$q75_pass <- is.finite(tt$q75) & tt$q75 > 0
  tt$threshold <- tt$q75_pass & !is.na(tt$qvalue) & tt$qvalue < 0.05 &
    is.finite(tt$logFC) & abs(tt$logFC) > 1 &
    is.finite(tt$max_mean_CPM) & tt$max_mean_CPM >= min_cpm
  tt[order(tt$qvalue, na.last = TRUE), , drop = FALSE]
}

make_unique_dash <- function(x) {
  x <- as.character(x)
  idx_within <- ave(seq_along(x), x, FUN = seq_along)
  dup <- duplicated(x) | duplicated(x, fromLast = TRUE)
  x[dup] <- paste0(x[dup], "-", idx_within[dup])
  x
}
stopifnot(!anyDuplicated(gene_anno$ensembl_gene_id))
miss <- is.na(gene_anno$hgnc_symbol) | gene_anno$hgnc_symbol == ""
gene_anno$hgnc_symbol[miss] <- gene_anno$ensembl_gene_id[miss]
gene_anno$hgnc_symbol_unique <- make_unique_dash(gene_anno$hgnc_symbol)

annotate_de_table <- function(tt) {
  id <- rownames(tt)
  i <- match(id, gene_anno$ensembl_gene_id)
  sym <- gene_anno$hgnc_symbol[i]
  symu <- gene_anno$hgnc_symbol_unique[i]
  sym[is.na(sym) | sym == ""] <- id[is.na(sym) | sym == ""]
  symu[is.na(symu) | symu == ""] <- id[is.na(symu) | symu == ""]
  cbind(data.frame(ensembl_gene_id = id, gene_symbol = sym, gene_symbol_unique = symu,
                   check.names = FALSE), tt)
}

# Load counts, align metadata, restrict to annotated genes and the q75 > 0
# universe (computed over ALL samples in the count matrix, as in the notebooks).
load_dataset <- function(count_dir, meta, id_col, strip_lane_suffix = TRUE) {
  cts <- fread(file.path(dracc_dir, count_dir, "Tables/genesCounts.csv"), data.table = FALSE)
  # GCM1 metadata file_id keeps the _S##_L### suffix; the other datasets drop it.
  if (strip_lane_suffix) colnames(cts) <- gsub("_S[0-9]+_L[0-9]+$", "", colnames(cts))
  stopifnot(setequal(meta[[id_col]], names(cts)[-1]), !anyDuplicated(meta[[id_col]]))
  meta <- meta[match(names(cts)[-1], meta[[id_col]]), ]
  cts_dat <- data.matrix(cts[, -1])
  rownames(cts_dat) <- cts$Name
  cts_dat <- cts_dat[rownames(cts_dat) %in% gene_anno$ensembl_gene_id, ]
  gene_info <- data.frame(q75 = apply(cts_dat, 1, quantile, probs = 0.75))
  keep <- gene_info$q75 > 0
  list(cts = cts_dat[keep, ], meta = meta, gene_info = gene_info[keep, , drop = FALSE])
}

# Fit one KO-vs-WT contrast. `clone` labels the KO clone of each sample.
run_contrast <- function(ds, samples, is_wt, clone, label) {
  counts <- ds$cts[, samples, drop = FALSE]
  group <- factor(ifelse(is_wt, "WT", "KO"), levels = c("WT", "KO"))
  block <- as.character(clone)
  if (wt_block_mode == "independent") {
    block[is_wt] <- paste0("WT_", seq_len(sum(is_wt)))
  } else {
    block[is_wt] <- "WT_parental"
  }
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("WT", "KO")
  y <- suppressMessages(calcNormFactors(DGEList(counts)))
  fit <- fit_limma_with_optional_duplicate_correlation(y, design, block = block, fit_label = label)
  cm <- makeContrasts(KOvsWT = KO - WT, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tt <- topTable(fit2, coef = "KOvsWT", number = Inf)
  tt <- add_qvalue_and_cpm(tt, counts, group, ds$gene_info)
  tt <- annotate_de_table(tt)
  write.table(tt, file.path(out_dir, paste0(label, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  data.frame(
    contrast = label, wt_block = wt_block_mode,
    n_WT = sum(is_wt), n_KO = sum(!is_wt),
    blocks = paste(sprintf("%s(%d)", names(table(block)), table(block)), collapse = "; "),
    dupcor_first_pass = fit$duplicate_correlation_first_pass,
    dupcor = fit$duplicate_correlation,
    dupcor_used = fit$duplicate_correlation_used,
    df_prior = fit2$df.prior[1],
    n_DE = sum(tt$threshold), n_up = sum(tt$threshold & tt$logFC > 0),
    n_down = sum(tt$threshold & tt$logFC < 0)
  )
}

summary_rows <- list()

# ---------------------------------------------------------------------------
# EOMES (MSK-processed data; definitive endoderm)
# ---------------------------------------------------------------------------
message("EOMES")
meta <- fread(file.path(meta_dir, "MSK_meta_data.tsv"), data.table = FALSE)
ds <- load_dataset("MSK_exchange_experiment_DRACC_processed_March_2026", meta, "Name")
for (ctr in list(list("JAX", c("JAX_KO_10", "JAX_KO_12")), list("MSK", c("MSK_KO_9", "MSK_KO_2")))) {
  keep <- ds$meta$Dcondition %in% c(ctr[[2]], "KOLF")
  m <- ds$meta[keep, ]
  summary_rows[[length(summary_rows) + 1]] <- run_contrast(
    ds, which(keep), m$ko_gene == "WT", m$Dcondition, paste0("EOMES_DE_", ctr[[1]]))
}

# ---------------------------------------------------------------------------
# GCM1 (JAX-processed data; trophoblast / syncytium)
# ---------------------------------------------------------------------------
message("GCM1")
meta <- fread(file.path(meta_dir, "JAX_meta_data.tsv"), data.table = FALSE)
meta$model_condition <- substr(meta$model_system, 1, nchar(meta$model_system) - 1)
ds <- load_dataset("JAX_exchange_experiment_DRACC_processed_August2025", meta, "file_id",
                   strip_lane_suffix = FALSE)
for (ctr in list(list("JAX", c("JAX_C01_Rep", "JAX_E04_Rep")), list("MSK", c("MSK_Clone1_Rep", "MSK_Clone2_Rep")))) {
  keep <- ds$meta$model_condition %in% c(ctr[[2]], "JAX_MSK_WT")
  m <- ds$meta[keep, ]
  is_wt <- m$ko_strategy == "CTRL_KO_WT" | m$ko_gene == "WT" | m$model_condition == "JAX_MSK_WT"
  clone <- sapply(strsplit(m$model_condition, "_"), `[`, 2)
  summary_rows[[length(summary_rows) + 1]] <- run_contrast(
    ds, which(keep), is_wt, clone, paste0("GCM1_DE_", ctr[[1]]))
}

# ---------------------------------------------------------------------------
# ISL1 (NW-processed data)
# ---------------------------------------------------------------------------
message("ISL1")
meta <- fread(file.path(meta_dir, "NW_meta_data.tsv"), data.table = FALSE)
ds <- load_dataset("NW_exchange_experiment_DRACC_processed_March2026", meta, "Name")
for (ctr in list(list("JAX", c("JAX_CLONE_AD2", "JAX_CLONE_DO4")), list("MSK", c("MSK_CLONE_05", "MSK_CLONE_07")))) {
  keep <- ds$meta$Dcondition %in% c(ctr[[2]], "KOLF_PARENT_NA")
  m <- ds$meta[keep, ]
  summary_rows[[length(summary_rows) + 1]] <- run_contrast(
    ds, which(keep), m$ko_gene == "WT", m$Dcondition, paste0("ISL1_DE_", ctr[[1]]))
}

# ---------------------------------------------------------------------------
# NKX2-1 (UCSF-processed data). MSK contrast uses only clone K03, as in the
# primary notebook.
# ---------------------------------------------------------------------------
message("NKX2-1")
meta <- fread(file.path(meta_dir, "UCSF_meta_data.tsv"), data.table = FALSE)
meta$KO_WT <- "WT"
meta$KO_WT[meta$Condition %in% c("full knockdown", "degron plus IAA", "CRISPRi plus guide plus Dox")] <- "KO"
ds <- load_dataset("UCSF_exchange_experiment_DRACC_processed_August205", meta, "Name")
for (ctr in list(list("JAX", c("JAX_Clone_A03", "JAX_Clone_A06")), list("MSK", c("MSK_Clone_K03")))) {
  keep <- ds$meta$Dcondition %in% c(ctr[[2]], "KOLF_Parent")
  m <- ds$meta[keep, ]
  summary_rows[[length(summary_rows) + 1]] <- run_contrast(
    ds, which(keep), m$KO_WT == "WT", m$Dcondition, paste0("NKX2-1_DE_", ctr[[1]]))
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file.path(out_dir, "model_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_df, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
