#!/usr/bin/env Rscript

# 2026.09.11 --------------------------------------------------------------
# TCGA-BRCA: TNBC vs. Non-TNBC biomarker discovery
# Human-readable analytic example for an agent skill.
#
# 1. QC before batch correction
# 2. TSS/SVA/RUVr/RUVs estimation, combination comparison and batch correction
# 3. DESeq2 differential-expression test
# 4. Volcano plot, DEG heatmap and gene-symbol tables
# 5. Hallmark and GO Biological Process GSEA
# 6. Gene-wise overall-survival analysis
# 7. Merge DE and survival results and write an HTML report
#
# This file deliberately shows the analysis as sequential R code, following
# 20250402_tcga_stad_analysis_manual_annotation_v4_hr.R. It does not source
# HR_biomarker_discovery/R/*.R. Edit the setting block below for a new study.


# 0. Packages and analysis settings --------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(sva)
  library(RUVSeq)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(clusterProfiler)
  library(survival)
  library(data.table)
  library(jsonlite)
})

## 0.1. Project and input paths ####
command_file <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (length(command_file)) {
  script_file <- normalizePath(sub('^--file=', '', command_file[[1]]), mustWork = TRUE)
  project_dir <- dirname(dirname(script_file))
} else {
  script_file <- NA_character_
  project_dir <- normalizePath('HR_biomarker_discovery', mustWork = TRUE)
}

default_cache_dir <- file.path(
  project_dir, '..', 'omics_pipeline', 'output', 'omicraft_core',
  '20260825_tcga_brca_tnbc_vs_non_tnbc'
)

counts_file <- Sys.getenv(
  'HR_COUNTS_RDS',
  file.path(default_cache_dir, '03_comparison', '20260825_counts_selected_raw_hr.rds')
)
metadata_file <- Sys.getenv(
  'HR_METADATA_RDS',
  file.path(default_cache_dir, '03_comparison', '20260825_metadata_selected_hr.rds')
)
annotation_file <- Sys.getenv(
  'HR_ANNOTATION_RDS',
  file.path(default_cache_dir, '03_comparison', '20260825_gene_annotation_selected_hr.rds')
)
clinical_file <- Sys.getenv(
  'HR_CLINICAL_RDS',
  file.path(default_cache_dir, '01_download', '20260825_BRCA_receptor_clinical_raw_hr.rds')
)
gene_sets_file <- Sys.getenv(
  'HR_GENE_SETS_RDS',
  file.path(default_cache_dir, '08_gsea', '20260825_MSigDB_gene_sets_snapshot_hr.rds')
)

run_id <- format(Sys.time(), '%Y%m%d_%H%M%S')
output_dir <- Sys.getenv(
  'HR_OUTPUT_DIR',
  file.path(project_dir, 'output', paste0(run_id, '_tcga_brca_tnbc_analytic'))
)

input_files <- c(
  counts = counts_file,
  metadata = metadata_file,
  annotation = annotation_file,
  clinical = clinical_file,
  gene_sets = gene_sets_file
)
if (any(!file.exists(input_files))) {
  stop('Missing input file(s): ', paste(input_files[!file.exists(input_files)], collapse = ', '))
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop('Output directory is not empty: ', output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

dir_qc <- file.path(output_dir, '01_qc')
dir_batch <- file.path(output_dir, '02_batch')
dir_de <- file.path(output_dir, '03_de')
dir_de_plots <- file.path(output_dir, '04_de_plots')
dir_gsea <- file.path(output_dir, '05_gsea')
dir_survival <- file.path(output_dir, '06_survival')
dir_merge <- file.path(output_dir, '07_merge')
for (directory in c(dir_qc, dir_batch, dir_de, dir_de_plots, dir_gsea, dir_survival, dir_merge)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

## 0.2. Analysis parameters ####
seed <- 20260911L
min_count <- 10L
min_samples <- 5L
pca_top_n <- 3000L
batch_estimation_genes <- 5000L
max_sva <- as.integer(Sys.getenv('HR_MAX_SVA', '2'))
max_ruvr <- as.integer(Sys.getenv('HR_MAX_RUVR', '2'))
max_ruvs <- as.integer(Sys.getenv('HR_MAX_RUVS', '2'))
max_latent_factors <- 3L
min_residual_df <- 10L
de_alpha <- 0.05
lfc_cutoff <- 1
heatmap_top_n <- 50L
volcano_label_n <- 12L
gsea_min_size <- 15L
gsea_max_size <- 500L
gsea_alpha <- 0.05
survival_strata <- c('pooled', 'TNBC', 'Non_TNBC')
save_large_objects <- tolower(Sys.getenv('HR_SAVE_LARGE_OBJECTS', 'true')) %in% c('1', 'true', 'yes')
set.seed(seed)

analysis_parameters <- list(
  seed = seed,
  comparison = 'TNBC vs Non_TNBC',
  positive_log2FC = 'Higher in TNBC',
  expression_filter = list(min_count = min_count, min_samples = min_samples),
  qc = list(pca_top_n = pca_top_n),
  batch = list(
    max_sva = max_sva, max_ruvr = max_ruvr, max_ruvs = max_ruvs,
    max_latent_factors = max_latent_factors,
    RUV_control_genes = 'all factor-estimation genes; exploratory primary-reference assumption',
    RUVs_replicates = 'within-condition biological samples; exploratory, not technical replicates'
  ),
  de = list(alpha = de_alpha, absolute_log2FC = lfc_cutoff),
  gsea = list(collections = c('Hallmark', 'GOBP'), min_size = gsea_min_size,
              max_size = gsea_max_size, alpha = gsea_alpha),
  survival = list(strata = survival_strata, split = 'HIGH >= median; LOW < median',
                  gene_universe = 'same expression-filtered genes as DESeq2')
)
write_json(analysis_parameters, file.path(output_dir, 'analysis_parameters.json'),
           pretty = TRUE, auto_unbox = TRUE)

personal_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(colour = 'black'),
      plot.title = element_text(hjust = 0.5, face = 'bold')
    )
}
group_colors <- c(TNBC = '#C43C39', Non_TNBC = '#377EB8')


# 0. Input preparation ----------------------------------------------------
cat("input_preparation\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 0.1. Read raw counts and metadata ####
cnts <- readRDS(counts_file)
metadata <- as.data.frame(readRDS(metadata_file))
gene_metadata <- as.data.frame(readRDS(annotation_file))
clinical <- as.data.frame(readRDS(clinical_file))

if (!'sample_id' %in% names(metadata)) metadata$sample_id <- metadata$sample_barcode
if (!'patient_id' %in% names(metadata)) metadata$patient_id <- metadata$patient_barcode
if (!'tss' %in% names(metadata)) metadata$tss <- metadata$tss_code
if (!'gene_symbol' %in% names(gene_metadata)) gene_metadata$gene_symbol <- gene_metadata$gene_name

metadata <- metadata[match(colnames(cnts), metadata$sample_id), , drop = FALSE]
gene_metadata <- gene_metadata[match(rownames(cnts), gene_metadata$gene_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
rownames(gene_metadata) <- gene_metadata$gene_id

metadata$group <- factor(metadata$group, levels = c('Non_TNBC', 'TNBC'))
gene_metadata$gene_symbol <- trimws(as.character(gene_metadata$gene_symbol))
gene_metadata$gene_symbol[is.na(gene_metadata$gene_symbol) | gene_metadata$gene_symbol == ''] <- NA_character_
gene_metadata$gene_label <- ifelse(is.na(gene_metadata$gene_symbol), gene_metadata$gene_id, gene_metadata$gene_symbol)

stopifnot(
  is.matrix(cnts), is.numeric(cnts), all(is.finite(cnts)), all(cnts >= 0),
  all(cnts == round(cnts)), !anyDuplicated(rownames(cnts)), !anyDuplicated(colnames(cnts)),
  !anyDuplicated(metadata$sample_id), !anyDuplicated(metadata$patient_id),
  !anyDuplicated(gene_metadata$gene_id),
  identical(colnames(cnts), metadata$sample_id),
  identical(rownames(cnts), gene_metadata$gene_id),
  all(c('Non_TNBC', 'TNBC') %in% metadata$group)
)

## 0.2. Reconstruct overall survival ####
if (all(c('bcr_patient_barcode', 'vital_status', 'death_days_to', 'last_contact_days_to') %in% names(clinical))) {
  clinical_match <- clinical[match(metadata$patient_id, clinical$bcr_patient_barcode), , drop = FALSE]
  vital_status <- tolower(trimws(clinical_match$vital_status))
  metadata$os_event <- ifelse(vital_status == 'dead', 1L,
                              ifelse(vital_status == 'alive', 0L, NA_integer_))
  death_days <- suppressWarnings(as.numeric(as.character(clinical_match$death_days_to)))
  contact_days <- suppressWarnings(as.numeric(as.character(clinical_match$last_contact_days_to)))
  metadata$os_time_days <- ifelse(metadata$os_event == 1L, death_days,
                                  ifelse(metadata$os_event == 0L, contact_days, NA_real_))
} else if (all(c('patient_id', 'os_event', 'os_time_days') %in% names(clinical))) {
  clinical_match <- clinical[match(metadata$patient_id, clinical$patient_id), , drop = FALSE]
  metadata$os_event <- suppressWarnings(as.numeric(as.character(clinical_match$os_event)))
  metadata$os_time_days <- suppressWarnings(as.numeric(as.character(clinical_match$os_time_days)))
} else {
  stop('Unsupported clinical schema')
}

metadata$os_eligible <- is.finite(metadata$os_time_days) & metadata$os_time_days > 0 &
  !is.na(metadata$os_event) & metadata$os_event %in% c(0, 1)

## 0.3. Use one expression filter for DE and survival ####
keep <- rowSums(cnts >= min_count) >= min_samples
if (sum(keep) < 50L) stop('Fewer than 50 genes pass the expression filter')
gene_metadata$de_filter_pass <- keep
cnts_de <- cnts[keep, , drop = FALSE]
gene_metadata_de <- gene_metadata[keep, , drop = FALSE]

dds_input <- DESeqDataSetFromMatrix(
  countData = cnts_de,
  colData = metadata,
  design = ~ group
)
dds_input <- estimateSizeFactors(dds_input)
size_factors <- sizeFactors(dds_input)
vsd_before <- varianceStabilizingTransformation(dds_input, blind = FALSE)
vst_before <- assay(vsd_before)

dim(cnts)
dim(cnts_de)
table(metadata$group)
table(metadata$os_eligible, metadata$group)

fwrite(metadata, file.path(output_dir, 'sample_metadata.tsv'), sep = '\t', na = 'NA')
fwrite(gene_metadata, file.path(output_dir, 'gene_annotation_and_filter.tsv'), sep = '\t', na = 'NA')
fwrite(data.frame(sample_id = names(size_factors), size_factor = size_factors),
       file.path(output_dir, 'size_factors.tsv'), sep = '\t')
saveRDS(cnts, file.path(output_dir, 'counts_all.rds'), compress = FALSE)
saveRDS(vst_before, file.path(output_dir, 'vst_before.rds'), compress = FALSE)

input_manifest <- list(
  input_files = as.list(normalizePath(input_files)),
  input_md5 = as.list(tools::md5sum(input_files)),
  samples = ncol(cnts), genes_raw = nrow(cnts), genes_de_and_survival = sum(keep),
  group_counts = as.list(table(metadata$group)),
  os_eligible = sum(metadata$os_eligible), os_events = sum(metadata$os_event[metadata$os_eligible])
)
write_json(input_manifest, file.path(output_dir, 'input_manifest.json'),
           pretty = TRUE, auto_unbox = TRUE, na = 'null')


# 1. QC before batch correction ------------------------------------------
cat("qc\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 1.1. PCA plot ####
gene_variance <- matrixStats::rowVars(vst_before)
pca_genes <- rownames(vst_before)[
  head(order(gene_variance, decreasing = TRUE), min(pca_top_n, sum(gene_variance > 0)))
]
pca_before_fit <- prcomp(t(vst_before[pca_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_before_percent <- 100 * pca_before_fit$sdev^2 / sum(pca_before_fit$sdev^2)
pca_before <- cbind(
  metadata[, c('sample_id', 'patient_id', 'group', 'tss')],
  as.data.frame(pca_before_fit$x[, 1:5, drop = FALSE])
)
fwrite(pca_before, file.path(dir_qc, 'before_batch_pca_scores.tsv'), sep = '\t')
fwrite(data.frame(gene_id = pca_genes), file.path(dir_qc, 'before_batch_pca_genes.tsv'), sep = '\t')

p_pca_group_before <- ggplot(pca_before, aes(PC1, PC2, colour = group)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = group_colors) +
  personal_theme() +
  labs(title = 'Before batch correction', colour = 'Group',
       x = sprintf('PC1 (%.1f%%)', pca_before_percent[1]),
       y = sprintf('PC2 (%.1f%%)', pca_before_percent[2]))
ggsave(file.path(dir_qc, 'before_batch_pca_group.pdf'), p_pca_group_before, width = 7, height = 6)
ggsave(file.path(dir_qc, 'before_batch_pca_group.png'), p_pca_group_before, width = 7, height = 6, dpi = 180)

p_pca_tss_before <- ggplot(pca_before, aes(PC1, PC2, colour = tss)) +
  geom_point(size = 1.8, alpha = 0.75) +
  guides(colour = 'none') + personal_theme() +
  labs(title = 'Before batch correction: Tissue Source Site',
       x = sprintf('PC1 (%.1f%%)', pca_before_percent[1]),
       y = sprintf('PC2 (%.1f%%)', pca_before_percent[2]))
ggsave(file.path(dir_qc, 'before_batch_pca_tss.pdf'), p_pca_tss_before, width = 7, height = 6)
ggsave(file.path(dir_qc, 'before_batch_pca_tss.png'), p_pca_tss_before, width = 7, height = 6, dpi = 180)

## 1.2. Sample distance plot ####
sample_distance_before <- dist(t(vst_before))
sample_distance_matrix_before <- as.matrix(sample_distance_before)
sample_annotation <- data.frame(group = metadata$group, row.names = metadata$sample_id)
sample_hclust_before <- hclust(sample_distance_before, method = 'ward.D')
saveRDS(sample_distance_matrix_before, file.path(dir_qc, 'before_batch_sample_distance.rds'))

pdf(file.path(dir_qc, 'before_batch_sample_distance.pdf'), width = 9, height = 8)
pheatmap(sample_distance_matrix_before, cluster_rows = sample_hclust_before,
         cluster_cols = sample_hclust_before, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = sample_annotation, annotation_colors = list(group = group_colors),
         border_color = NA, main = 'Before batch correction: sample distance')
dev.off()
png(file.path(dir_qc, 'before_batch_sample_distance.png'), width = 1600, height = 1450, res = 180)
pheatmap(sample_distance_matrix_before, cluster_rows = sample_hclust_before,
         cluster_cols = sample_hclust_before, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = sample_annotation, annotation_colors = list(group = group_colors),
         border_color = NA, main = 'Before batch correction: sample distance')
dev.off()


# 2. Batch correction -----------------------------------------------------
cat("batch_correction\n", file=file.path(output_dir, "progress.log"), append=TRUE)

group_design <- model.matrix(~ group, metadata)
batch_factors <- data.frame(row.names = metadata$sample_id)
factor_status <- list()

## 2.1. Tissue Source Site ####
raw_tss <- as.character(metadata$tss)
raw_tss[is.na(raw_tss) | raw_tss == ''] <- 'TSS_UNKNOWN'
rare_tss <- names(which(table(raw_tss) < 3L))
model_tss <- raw_tss
model_tss[model_tss %in% setdiff(rare_tss, 'TSS_UNKNOWN')] <- 'TSS_OTHER'
batch_factors$TSS <- droplevels(factor(model_tss))
if (nlevels(batch_factors$TSS) < 2L) batch_factors$TSS <- NULL
factor_status[['TSS']] <- data.frame(
  method = 'TSS', status = if ('TSS' %in% names(batch_factors)) 'available' else 'excluded',
  n_factors = 0L, note = 'Sites with n<3 pooled as TSS_OTHER'
)
fwrite(data.frame(sample_id = metadata$sample_id, original_tss = raw_tss, model_tss = model_tss),
       file.path(dir_batch, 'tss_mapping.tsv'), sep = '\t')

## 2.2. Select variable genes for factor estimation ####
batch_gene_order <- order(matrixStats::rowVars(vst_before), decreasing = TRUE)
batch_gene_ids <- rownames(vst_before)[head(batch_gene_order, min(batch_estimation_genes, nrow(vst_before)))]
batch_counts <- cnts_de[batch_gene_ids, , drop = FALSE]
batch_norm_counts <- sweep(batch_counts, 2L, size_factors[colnames(batch_counts)], '/')
fwrite(data.frame(gene_id = batch_gene_ids), file.path(dir_batch, 'factor_estimation_genes.tsv'), sep = '\t')

## 2.3. SVA ####
sva_fit <- NULL
sva_error <- NA_character_
if (max_sva > 0L) {
  n_sva <- min(max_sva, ncol(batch_counts) - ncol(group_design) - min_residual_df)
  if (n_sva > 0L) {
    sva_fit <- tryCatch(
      svaseq(batch_norm_counts, mod = group_design, mod0 = model.matrix(~ 1, metadata),
             n.sv = n_sva, B = 5),
      error = function(error) {
        sva_error <<- conditionMessage(error)
        NULL
      }
    )
  } else {
    sva_error <- 'Insufficient residual degrees of freedom for SVA'
  }
  if (!is.null(sva_fit)) {
    sva_matrix <- scale(sva_fit$sv)
    colnames(sva_matrix) <- paste0('SV', seq_len(ncol(sva_matrix)))
    rownames(sva_matrix) <- metadata$sample_id
    batch_factors <- cbind(batch_factors, as.data.frame(sva_matrix))
  }
}
factor_status[['SVA']] <- data.frame(
  method = 'SVA', status = if (!is.null(sva_fit)) 'available' else if (max_sva == 0L) 'disabled' else 'excluded',
  n_factors = if (is.null(sva_fit)) 0L else ncol(sva_fit$sv),
  note = if (is.null(sva_fit)) ifelse(is.na(sva_error), paste0('max_sva=', max_sva), sva_error)
         else 'svaseq normalized counts; full ~group, null ~1'
)

## 2.4. RUVr ####
ruvr_fit <- NULL
ruvr_error <- NA_character_
ruv_residuals <- NULL
if (max_ruvr > 0L || max_ruvs > 0L) {
  ruv_residuals <- tryCatch({
    edge_object <- DGEList(counts = batch_counts, group = metadata$group)
    edge_object <- calcNormFactors(edge_object, method = 'upperquartile')
    edge_object <- estimateGLMCommonDisp(edge_object, group_design)
    edge_object <- estimateGLMTagwiseDisp(edge_object, group_design)
    edge_fit <- glmFit(edge_object, group_design)
    residuals(edge_fit, type = 'deviance')
  }, error = function(error) {
    ruvr_error <<- conditionMessage(error)
    NULL
  })
}
if (max_ruvr > 0L && !is.null(ruv_residuals)) {
  ruvr_fit <- tryCatch(
    RUVr(batch_counts, cIdx = seq_len(nrow(batch_counts)), k = max_ruvr,
         residuals = ruv_residuals),
    error = function(error) {
      ruvr_error <<- conditionMessage(error)
      NULL
    }
  )
  if (!is.null(ruvr_fit)) {
    ruvr_matrix <- scale(ruvr_fit$W)
    colnames(ruvr_matrix) <- paste0('RUVr', seq_len(ncol(ruvr_matrix)))
    rownames(ruvr_matrix) <- metadata$sample_id
    batch_factors <- cbind(batch_factors, as.data.frame(ruvr_matrix))
  }
}
factor_status[['RUVr']] <- data.frame(
  method = 'RUVr', status = if (!is.null(ruvr_fit)) 'available' else if (max_ruvr == 0L) 'disabled' else 'excluded',
  n_factors = if (is.null(ruvr_fit)) 0L else ncol(ruvr_fit$W),
  note = if (is.null(ruvr_fit)) ifelse(is.na(ruvr_error), paste0('max_ruvr=', max_ruvr), ruvr_error)
         else 'Upper-quartile edgeR group-model deviance residuals; all estimation genes as controls'
)

## 2.5. RUVs ####
# makeGroups(group) treats within-condition biological samples as the
# reference replicate structure. This follows the primary code but is an
# exploratory assumption; these samples are not technical replicates.
ruvs_fit <- NULL
ruvs_error <- NA_character_
if (max_ruvs > 0L) {
  ruvs_fit <- tryCatch(
    RUVs(batch_counts, cIdx = seq_len(nrow(batch_counts)), k = max_ruvs,
         scIdx = makeGroups(metadata$group)),
    error = function(error) {
      ruvs_error <<- conditionMessage(error)
      NULL
    }
  )
  if (!is.null(ruvs_fit)) {
    ruvs_matrix <- scale(ruvs_fit$W)
    colnames(ruvs_matrix) <- paste0('RUVs', seq_len(ncol(ruvs_matrix)))
    rownames(ruvs_matrix) <- metadata$sample_id
    batch_factors <- cbind(batch_factors, as.data.frame(ruvs_matrix))
  }
}
factor_status[['RUVs']] <- data.frame(
  method = 'RUVs', status = if (!is.null(ruvs_fit)) 'available' else if (max_ruvs == 0L) 'disabled' else 'excluded',
  n_factors = if (is.null(ruvs_fit)) 0L else ncol(ruvs_fit$W),
  note = if (is.null(ruvs_fit)) ifelse(is.na(ruvs_error), paste0('max_ruvs=', max_ruvs), ruvs_error)
         else 'Exploratory within-condition reference; not technical replicates'
)

factor_status_table <- do.call(rbind, factor_status)
rownames(factor_status_table) <- NULL
fwrite(factor_status_table, file.path(dir_batch, 'factor_status.tsv'), sep = '\t', na = 'NA')
fwrite(cbind(sample_id = rownames(batch_factors), batch_factors),
       file.path(dir_batch, 'all_batch_covariates.tsv'), sep = '\t')

## 2.6. Compare valid factor combinations ####
batch_model_data <- cbind(metadata, batch_factors)
batch_grid <- expand.grid(
  tss = 0:1, sva = 0:max_sva, ruvr = 0:max_ruvr, ruvs = 0:max_ruvs,
  KEEP.OUT.ATTRS = FALSE
)

batch_pca_reference <- prcomp(t(vst_before[pca_genes, , drop = FALSE]),
                              center = TRUE, scale. = FALSE, rank. = 10L)
technical_tss_matrix <- if ('TSS' %in% names(batch_factors)) {
  model.matrix(~ TSS, batch_model_data)[, -1, drop = FALSE]
} else NULL

candidate_rows <- vector('list', nrow(batch_grid))
candidate_columns <- vector('list', nrow(batch_grid))

for (candidate_index in seq_len(nrow(batch_grid))) {
  grid_row <- batch_grid[candidate_index, ]
  batch_columns <- c(
    if (grid_row$tss == 1L) 'TSS',
    if (grid_row$sva > 0L) paste0('SV', seq_len(grid_row$sva)),
    if (grid_row$ruvr > 0L) paste0('RUVr', seq_len(grid_row$ruvr)),
    if (grid_row$ruvs > 0L) paste0('RUVs', seq_len(grid_row$ruvs))
  )
  candidate_id <- if (length(batch_columns)) paste(batch_columns, collapse = '+') else 'none'
  candidate_columns[[candidate_index]] <- batch_columns
  latent_factors <- grid_row$sva + grid_row$ruvr + grid_row$ruvs
  missing_factors <- setdiff(batch_columns, names(batch_factors))

  candidate_result <- data.frame(
    candidate_id = candidate_id,
    design_formula = paste('~', paste(c(batch_columns, 'group'), collapse = ' + ')),
    eligible = FALSE, reason = '', nuisance_columns = NA_integer_, latent_factors = latent_factors,
    residual_df = NA_integer_, condition_number = NA_real_, group_nuisance_r2 = NA_real_,
    tss_partial_r2 = NA_real_, score = NA_real_, stringsAsFactors = FALSE
  )

  if (length(missing_factors)) {
    candidate_result$reason <- paste0('unavailable_factors:', paste(missing_factors, collapse = ','))
  } else if (latent_factors > max_latent_factors) {
    candidate_result$reason <- 'latent_factor_limit'
  } else {
    candidate_formula <- as.formula(candidate_result$design_formula)
    candidate_design <- model.matrix(candidate_formula, batch_model_data)
    nuisance_design <- if (length(batch_columns)) {
      model.matrix(reformulate(batch_columns), batch_model_data)
    } else {
      matrix(1, nrow(metadata), 1L, dimnames = list(metadata$sample_id, '(Intercept)'))
    }
    candidate_result$nuisance_columns <- ncol(nuisance_design) - 1L
    candidate_result$residual_df <- nrow(candidate_design) - qr(candidate_design)$rank
    candidate_result$condition_number <- kappa(candidate_design)

    group_numeric <- as.numeric(metadata$group) - 1
    group_sse <- sum(qr.resid(qr(nuisance_design), group_numeric)^2)
    group_total <- sum((group_numeric - mean(group_numeric))^2)
    candidate_result$group_nuisance_r2 <- max(0, 1 - group_sse / group_total)

    if (qr(candidate_design)$rank < ncol(candidate_design)) {
      candidate_result$reason <- 'rank_deficient_design'
    } else if (candidate_result$residual_df < min_residual_df) {
      candidate_result$reason <- 'insufficient_residual_df'
    } else if (!is.finite(candidate_result$condition_number) || candidate_result$condition_number > 1e6) {
      candidate_result$reason <- 'ill_conditioned_design'
    } else if (candidate_result$group_nuisance_r2 > 0.95) {
      candidate_result$reason <- 'group_nuisance_confounding_limit'
    } else {
      candidate_adjusted_vst <- if (length(batch_columns)) {
        removeBatchEffect(vst_before, covariates = nuisance_design[, -1, drop = FALSE],
                          design = group_design)
      } else {
        vst_before
      }

      if (!is.null(technical_tss_matrix)) {
        centered_adjusted <- sweep(t(candidate_adjusted_vst[pca_genes, , drop = FALSE]),
                                   2L, batch_pca_reference$center, '-')
        adjusted_scores <- centered_adjusted %*% batch_pca_reference$rotation
        sse_group <- sum(qr.resid(qr(group_design), adjusted_scores)^2)
        sse_group_tss <- sum(qr.resid(qr(cbind(group_design, technical_tss_matrix)), adjusted_scores)^2)
        candidate_result$tss_partial_r2 <- if (sse_group < 1e-12) 0 else {
          max(0, min(1, (sse_group - sse_group_tss) / sse_group))
        }
      }

      technical_score <- ifelse(is.finite(candidate_result$tss_partial_r2),
                                candidate_result$tss_partial_r2, 0)
      confounding_penalty <- -log(max(1e-12, 1 - candidate_result$group_nuisance_r2))
      candidate_result$score <- technical_score +
        0.001 * candidate_result$nuisance_columns + 0.001 * confounding_penalty
      candidate_result$eligible <- TRUE
    }
  }
  candidate_rows[[candidate_index]] <- candidate_result
}

batch_candidates <- do.call(rbind, candidate_rows)
rownames(batch_candidates) <- NULL
eligible_candidates <- which(batch_candidates$eligible)
if (!length(eligible_candidates)) stop('No batch combination passed the design checks')

selected_index <- eligible_candidates[
  order(batch_candidates$score[eligible_candidates],
        batch_candidates$nuisance_columns[eligible_candidates],
        batch_candidates$candidate_id[eligible_candidates])[1]
]
none_index <- which(batch_candidates$candidate_id == 'none' & batch_candidates$eligible)
if (length(none_index) == 1L &&
    batch_candidates$score[none_index] - batch_candidates$score[selected_index] < 0.005) {
  selected_index <- none_index
}

batch_candidates$selected <- seq_len(nrow(batch_candidates)) == selected_index
selected_batch_id <- batch_candidates$candidate_id[selected_index]
selected_batch_columns <- candidate_columns[[selected_index]]
selected_design_formula <- as.formula(batch_candidates$design_formula[selected_index])
selected_batch_covariates <- batch_factors[, selected_batch_columns, drop = FALSE]

selected_nuisance_design <- if (length(selected_batch_columns)) {
  model.matrix(reformulate(selected_batch_columns), batch_model_data)
} else {
  matrix(1, nrow(metadata), 1L)
}
vst_batch_corrected <- if (length(selected_batch_columns)) {
  removeBatchEffect(vst_before, covariates = selected_nuisance_design[, -1, drop = FALSE],
                    design = group_design)
} else {
  vst_before
}

fwrite(batch_candidates, file.path(dir_batch, 'batch_candidates.tsv'), sep = '\t', na = 'NA')
fwrite(cbind(sample_id = metadata$sample_id, selected_batch_covariates),
       file.path(dir_batch, 'selected_batch_covariates.tsv'), sep = '\t')
writeLines(c(
  paste0('Selected combination: ', selected_batch_id),
  paste0('DESeq2 design: ', paste(deparse(selected_design_formula), collapse = ' ')),
  'Selection score uses residual TSS association plus complexity/confounding penalties.',
  'DEG, GSEA and survival outcomes are not used to select the batch model.',
  'removeBatchEffect output is used only for visualization; DESeq2 uses raw counts.'
), file.path(dir_batch, 'batch_selection_notes.txt'))

batch_candidates[batch_candidates$selected, ]
factor_status_table

## 2.7. PCA after batch correction ####
pca_after_fit <- prcomp(t(vst_batch_corrected[pca_genes, , drop = FALSE]),
                        center = TRUE, scale. = FALSE)
pca_after_percent <- 100 * pca_after_fit$sdev^2 / sum(pca_after_fit$sdev^2)
pca_after <- cbind(
  metadata[, c('sample_id', 'patient_id', 'group', 'tss')],
  as.data.frame(pca_after_fit$x[, 1:5, drop = FALSE])
)
fwrite(pca_after, file.path(dir_batch, 'after_batch_pca_scores.tsv'), sep = '\t')

p_pca_group_after <- ggplot(pca_after, aes(PC1, PC2, colour = group)) +
  geom_point(size = 1.8, alpha = 0.75) + scale_colour_manual(values = group_colors) +
  personal_theme() + labs(title = paste('After batch correction:', selected_batch_id),
                          x = sprintf('PC1 (%.1f%%)', pca_after_percent[1]),
                          y = sprintf('PC2 (%.1f%%)', pca_after_percent[2]))
ggsave(file.path(dir_batch, 'after_batch_pca_group.pdf'), p_pca_group_after, width = 7, height = 6)
ggsave(file.path(dir_batch, 'after_batch_pca_group.png'), p_pca_group_after, width = 7, height = 6, dpi = 180)

p_pca_tss_after <- ggplot(pca_after, aes(PC1, PC2, colour = tss)) +
  geom_point(size = 1.8, alpha = 0.75) + guides(colour = 'none') + personal_theme() +
  labs(title = paste('After batch correction:', selected_batch_id),
       x = sprintf('PC1 (%.1f%%)', pca_after_percent[1]),
       y = sprintf('PC2 (%.1f%%)', pca_after_percent[2]))
ggsave(file.path(dir_batch, 'after_batch_pca_tss.pdf'), p_pca_tss_after, width = 7, height = 6)
ggsave(file.path(dir_batch, 'after_batch_pca_tss.png'), p_pca_tss_after, width = 7, height = 6, dpi = 180)

## 2.8. Sample distance after batch correction ####
sample_distance_after <- dist(t(vst_batch_corrected))
sample_distance_matrix_after <- as.matrix(sample_distance_after)
sample_hclust_after <- hclust(sample_distance_after, method = 'ward.D')
saveRDS(sample_distance_matrix_after, file.path(dir_batch, 'after_batch_sample_distance.rds'))

pdf(file.path(dir_batch, 'after_batch_sample_distance.pdf'), width = 9, height = 8)
pheatmap(sample_distance_matrix_after, cluster_rows = sample_hclust_after,
         cluster_cols = sample_hclust_after, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = sample_annotation, annotation_colors = list(group = group_colors),
         border_color = NA, main = paste('After batch correction:', selected_batch_id))
dev.off()
png(file.path(dir_batch, 'after_batch_sample_distance.png'), width = 1600, height = 1450, res = 180)
pheatmap(sample_distance_matrix_after, cluster_rows = sample_hclust_after,
         cluster_cols = sample_hclust_after, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = sample_annotation, annotation_colors = list(group = group_colors),
         border_color = NA, main = paste('After batch correction:', selected_batch_id))
dev.off()


# 3. DESeq2 ---------------------------------------------------------------
cat("differential_expression\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 3.1. Add selected covariates to colData ####
de_metadata <- metadata
if (length(selected_batch_columns)) {
  for (batch_column in selected_batch_columns) {
    de_metadata[[batch_column]] <- selected_batch_covariates[[batch_column]]
  }
}
de_design_matrix <- model.matrix(selected_design_formula, de_metadata)
if (qr(de_design_matrix)$rank < ncol(de_design_matrix)) stop('Selected DE design is rank deficient')

## 3.2. Run DESeq2 using raw integer counts ####
dds <- DESeqDataSetFromMatrix(
  countData = cnts_de,
  colData = de_metadata,
  design = selected_design_formula
)
sizeFactors(dds) <- size_factors[colnames(dds)]
dds <- DESeq(dds, fitType = 'parametric', parallel = FALSE)
resultsNames(dds)

res <- as.data.frame(results(
  dds,
  contrast = c('group', 'TNBC', 'Non_TNBC'),
  alpha = de_alpha
))
res$gene_id <- rownames(res)
res$beta_converged <- as.logical(S4Vectors::mcols(dds)$betaConv)
res$de_status <- ifelse(!res$beta_converged, 'not_converged',
                        ifelse(is.na(res$pvalue), 'outlier_or_untestable',
                               ifelse(is.na(res$padj), 'independent_filtered', 'tested')))

## 3.3. LFC shrinkage ####
res$log2FC_apeglm <- NA_real_
res$lfcSE_apeglm <- NA_real_
res_shrunken <- NULL
if (requireNamespace('apeglm', quietly = TRUE) && 'group_TNBC_vs_Non_TNBC' %in% resultsNames(dds)) {
  res_shrunken <- lfcShrink(dds, coef = 'group_TNBC_vs_Non_TNBC', type = 'apeglm')
  res$log2FC_apeglm <- res_shrunken$log2FoldChange
  res$lfcSE_apeglm <- res_shrunken$lfcSE
}

## 3.4. Annotate and classify DE results ####
res$direction <- 'Not_DEG'
is_deg <- res$de_status == 'tested' & !is.na(res$padj) & res$padj < de_alpha &
  abs(res$log2FoldChange) >= lfc_cutoff
res$direction[is_deg & res$log2FoldChange > 0] <- 'Up_in_TNBC'
res$direction[is_deg & res$log2FoldChange < 0] <- 'Down_in_TNBC'

de_table <- merge(gene_metadata, res, by = 'gene_id', all.x = TRUE, sort = FALSE)
de_table <- de_table[match(gene_metadata$gene_id, de_table$gene_id), , drop = FALSE]
de_table$de_status[!de_table$de_filter_pass] <- 'count_filtered'
de_table$direction[!de_table$de_filter_pass] <- 'Not_tested'
de_table <- de_table[, c('gene_symbol', 'gene_id', setdiff(names(de_table), c('gene_symbol', 'gene_id')))]

deg_table <- de_table[de_table$direction %in% c('Up_in_TNBC', 'Down_in_TNBC'), , drop = FALSE]
deg_table <- deg_table[order(deg_table$padj, -abs(deg_table$log2FoldChange)), , drop = FALSE]

fwrite(de_table, file.path(dir_de, 'DE_all_genes.tsv'), sep = '\t', na = 'NA')
fwrite(deg_table, file.path(dir_de, 'DEG.tsv'), sep = '\t', na = 'NA')
fwrite(de_metadata, file.path(dir_de, 'design_metadata.tsv'), sep = '\t', na = 'NA')
fwrite(cbind(sample_id = rownames(de_design_matrix), as.data.frame(de_design_matrix)),
       file.path(dir_de, 'design_matrix.tsv'), sep = '\t')
saveRDS(de_table, file.path(dir_de, 'de_table.rds'))
if (save_large_objects) saveRDS(dds, file.path(dir_de, 'dds_fitted.rds'))

de_summary <- list(
  contrast = 'TNBC vs Non_TNBC', positive_log2FC = 'Higher in TNBC',
  formula = paste(deparse(selected_design_formula), collapse = ' '),
  alpha = de_alpha, absolute_log2FC = lfc_cutoff,
  direction_counts = as.list(table(de_table$direction)),
  status_counts = as.list(table(de_table$de_status))
)
write_json(de_summary, file.path(dir_de, 'de_summary.json'), pretty = TRUE, auto_unbox = TRUE)
table(de_table$direction)
head(deg_table[, c('gene_symbol', 'gene_id', 'baseMean', 'log2FoldChange', 'padj', 'direction')])


# 4. DE result plots ------------------------------------------------------
cat("de_visualizations\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 4.1. Volcano plot ####
volcano_data <- de_table[
  de_table$de_filter_pass & de_table$beta_converged %in% TRUE &
    is.finite(de_table$log2FoldChange) & !is.na(de_table$padj),
  , drop = FALSE
]
volcano_data$minus_log10_padj <- -log10(pmax(volcano_data$padj, .Machine$double.xmin))
volcano_data$plot_y <- pmin(volcano_data$minus_log10_padj, 100)
volcano_labels <- volcano_data[volcano_data$direction != 'Not_DEG', , drop = FALSE]
volcano_labels <- head(volcano_labels[order(volcano_labels$padj), ], volcano_label_n)

volcano_plot <- ggplot(volcano_data, aes(log2FoldChange, plot_y, colour = direction)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
  geom_hline(yintercept = -log10(de_alpha), linetype = 2) +
  scale_colour_manual(values = c(Up_in_TNBC = '#C43C39', Down_in_TNBC = '#54AEE7', Not_DEG = 'grey70')) +
  personal_theme() +
  labs(title = 'TNBC vs. Non-TNBC', x = 'log2 fold change (TNBC / Non-TNBC)',
       y = '-log10 adjusted p-value', colour = NULL,
       subtitle = paste(sum(volcano_data$direction == 'Up_in_TNBC'), 'up /',
                        sum(volcano_data$direction == 'Down_in_TNBC'), 'down'))
if (nrow(volcano_labels)) {
  volcano_plot <- volcano_plot +
    geom_text_repel(data = volcano_labels, aes(label = gene_label),
                    size = 3, max.overlaps = Inf, seed = seed)
}
ggsave(file.path(dir_de_plots, 'volcano.pdf'), volcano_plot, width = 8, height = 7)
ggsave(file.path(dir_de_plots, 'volcano.png'), volcano_plot, width = 8, height = 7, dpi = 180)
fwrite(volcano_data, file.path(dir_de_plots, 'volcano_data.tsv'), sep = '\t', na = 'NA')

## 4.2. DEG heatmap ####
heatmap_genes <- deg_table[!duplicated(deg_table$gene_label), , drop = FALSE]
heatmap_genes <- head(heatmap_genes, heatmap_top_n)
heatmap_ids <- intersect(heatmap_genes$gene_id, rownames(vst_batch_corrected))

if (length(heatmap_ids)) {
  heatmap_matrix <- vst_batch_corrected[heatmap_ids, , drop = FALSE]
  rownames(heatmap_matrix) <- heatmap_genes$gene_label[match(heatmap_ids, heatmap_genes$gene_id)]
  fwrite(cbind(gene_label = rownames(heatmap_matrix), as.data.frame(heatmap_matrix)),
         file.path(dir_de_plots, 'DEG_heatmap_expression.tsv'), sep = '\t')
  fwrite(heatmap_genes[match(heatmap_ids, heatmap_genes$gene_id), ],
         file.path(dir_de_plots, 'DEG_heatmap_selected_genes.tsv'), sep = '\t', na = 'NA')

  pdf(file.path(dir_de_plots, 'DEG_heatmap.pdf'), width = 12, height = 10)
  pheatmap(heatmap_matrix, scale = 'row', show_colnames = FALSE,
           annotation_col = sample_annotation, annotation_colors = list(group = group_colors),
           border_color = NA, fontsize_row = 7,
           main = paste('Top', length(heatmap_ids), 'DEGs'))
  dev.off()
  png(file.path(dir_de_plots, 'DEG_heatmap.png'), width = 2160, height = 1800, res = 180)
  pheatmap(heatmap_matrix, scale = 'row', show_colnames = FALSE,
           annotation_col = sample_annotation, annotation_colors = list(group = group_colors),
           border_color = NA, fontsize_row = 7,
           main = paste('Top', length(heatmap_ids), 'DEGs'))
  dev.off()
} else {
  write_json(list(status = 'not_available', reason = 'No genes meet the DEG threshold'),
             file.path(dir_de_plots, 'DEG_heatmap_status.json'), auto_unbox = TRUE)
  fwrite(data.frame(gene_id = character(), gene_symbol = character()),
         file.path(dir_de_plots, 'DEG_heatmap_selected_genes.tsv'), sep = '\t')
}

rm(dds, dds_input, vsd_before, res_shrunken, pca_before_fit, pca_after_fit,
   sample_distance_matrix_before, sample_distance_matrix_after,
   batch_counts, batch_norm_counts, ruv_residuals)
gc(verbose = FALSE)


# 5. GSEA -----------------------------------------------------------------
cat("gsea\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 5.1. Rank all valid DE-tested genes ####
gsea_input <- de_table[
  !is.na(de_table$gene_symbol) & de_table$gene_symbol != '' &
    is.finite(de_table$stat) & de_table$beta_converged %in% TRUE,
  c('gene_id', 'gene_symbol', 'stat')
]
gsea_ranking_table <- aggregate(gsea_input$stat,
                                list(gene_symbol = gsea_input$gene_symbol), mean)
names(gsea_ranking_table)[2] <- 'mean_Wald_stat'
gsea_ranking_table <- gsea_ranking_table[
  order(-gsea_ranking_table$mean_Wald_stat, gsea_ranking_table$gene_symbol),
]
ranks <- setNames(gsea_ranking_table$mean_Wald_stat, gsea_ranking_table$gene_symbol)
fwrite(gsea_ranking_table, file.path(dir_gsea, 'GSEA_symbol_ranking.tsv'), sep = '\t')
fwrite(gsea_input, file.path(dir_gsea, 'GSEA_gene_id_to_symbol_audit.tsv'), sep = '\t')

## 5.2. Hallmark and GO Biological Process ####
gene_sets <- readRDS(gene_sets_file)
if (!all(c('Hallmark', 'GOBP') %in% names(gene_sets))) {
  stop('Gene-set snapshot must contain Hallmark and GOBP')
}
BiocParallel::register(BiocParallel::SerialParam())
gsea_summary <- list()

for (collection in c('Hallmark', 'GOBP')) {
  collection_data <- gene_sets[[collection]]
  term2gene <- unique(as.data.frame(collection_data$term2gene)[, 1:2])
  names(term2gene) <- c('term', 'gene')
  term2name <- unique(as.data.frame(collection_data$term2name)[, 1:2])
  names(term2name) <- c('term', 'name')

  set.seed(seed)
  gsea_object <- GSEA(
    geneList = ranks,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    minGSSize = gsea_min_size,
    maxGSSize = gsea_max_size,
    pvalueCutoff = 1,
    pAdjustMethod = 'BH',
    eps = 0,
    seed = TRUE,
    verbose = FALSE,
    by = 'fgsea'
  )
  gsea_table <- as.data.frame(gsea_object)
  gsea_significant <- gsea_table[
    !is.na(gsea_table$p.adjust) & gsea_table$p.adjust < gsea_alpha,
    , drop = FALSE
  ]
  fwrite(gsea_table, file.path(dir_gsea, paste0(collection, '_GSEA_all.tsv')),
         sep = '\t', na = 'NA')
  fwrite(gsea_significant, file.path(dir_gsea, paste0(collection, '_GSEA_significant.tsv')),
         sep = '\t', na = 'NA')
  saveRDS(gsea_object, file.path(dir_gsea, paste0(collection, '_GSEA.rds')))

  gsea_significant <- gsea_significant[
    order(gsea_significant$p.adjust, -abs(gsea_significant$NES)),
    , drop = FALSE
  ]
  gsea_top <- rbind(
    head(gsea_significant[gsea_significant$NES > 0, , drop = FALSE], 10L),
    head(gsea_significant[gsea_significant$NES < 0, , drop = FALSE], 10L)
  )

  if (nrow(gsea_top)) {
    gsea_top$pathway <- factor(make.unique(gsub('_', ' ', gsea_top$Description)),
                               levels = rev(make.unique(gsub('_', ' ', gsea_top$Description))))
    gsea_top$direction <- ifelse(gsea_top$NES > 0, 'TNBC', 'Non_TNBC')
    nes_plot <- ggplot(gsea_top, aes(NES, pathway, fill = direction)) +
      geom_col() + geom_vline(xintercept = 0) +
      scale_fill_manual(values = group_colors) + personal_theme() +
      labs(title = paste(collection, 'significant pathways'), y = NULL, fill = 'Enriched in')
  } else {
    write_json(list(status = 'not_available', reason = 'No significant pathways'),
               file.path(dir_gsea, paste0(collection, '_NES_plot_status.json')), auto_unbox = TRUE)
  }
  if (exists('gsea_top') && nrow(gsea_top) > 0) {
  ggsave(file.path(dir_gsea, paste0(collection, '_NES_top.pdf')), nes_plot,
         width = 13, height = max(5, 2 + 0.3 * nrow(gsea_top)))
  ggsave(file.path(dir_gsea, paste0(collection, '_NES_top.png')), nes_plot,
         width = 13, height = max(5, 2 + 0.3 * nrow(gsea_top)), dpi = 180)
  }

  gsea_summary[[collection]] <- list(
    tested = nrow(gsea_table), significant = nrow(gsea_significant),
    database = collection_data$provenance
  )
}

write_json(list(
  rank_metric = 'Mean DESeq2 Wald statistic per gene symbol; no DEG filter',
  positive_NES = 'TNBC', ranked_symbols = length(ranks),
  collections = gsea_summary
), file.path(dir_gsea, 'gsea_summary.json'), pretty = TRUE, auto_unbox = TRUE, na = 'null')


# 6. Survival analysis ----------------------------------------------------
cat("survival\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 6.1. Prepare the exact DE-input gene set ####
norm.cts <- sweep(cnts_de, 2L, size_factors[colnames(cnts_de)], '/')
stopifnot(
  identical(rownames(norm.cts), gene_metadata_de$gene_id),
  identical(rownames(norm.cts), res$gene_id),
  identical(colnames(norm.cts), metadata$sample_id)
)

valid_survival <- metadata$os_eligible & !is.na(metadata$group)
survival_rows <- vector('list', length(survival_strata) * nrow(norm.cts))
survival_row_index <- 0L

## 6.2. Median split, log-rank test and Cox HR ####
for (stratum in survival_strata) {
  in_stratum <- if (stratum == 'pooled') {
    rep(TRUE, nrow(metadata))
  } else {
    metadata$group == stratum
  }
  in_stratum[is.na(in_stratum)] <- FALSE
  analysis_samples <- valid_survival & in_stratum

  for (gene_index in seq_len(nrow(norm.cts))) {
    survival_row_index <- survival_row_index + 1L
    gene_expression_all <- norm.cts[gene_index, in_stratum]
    median_expression <- median(gene_expression_all[is.finite(gene_expression_all)], na.rm = TRUE)
    gene_expression <- norm.cts[gene_index, analysis_samples]
    survival_group <- factor(ifelse(gene_expression >= median_expression, 'HIGH', 'LOW'),
                             levels = c('LOW', 'HIGH'))
    survival_data <- data.frame(
      time = metadata$os_time_days[analysis_samples],
      event = metadata$os_event[analysis_samples],
      expression_group = survival_group
    )

    n_low <- sum(survival_group == 'LOW', na.rm = TRUE)
    n_high <- sum(survival_group == 'HIGH', na.rm = TRUE)
    events <- sum(survival_data$event)
    logrank_chisq <- logrank_p <- hr <- hr_lower95 <- hr_upper95 <- cox_p <- NA_real_
    status <- 'ok'
    warning_text <- ''

    if (!is.finite(median_expression)) {
      status <- 'no_finite_expression'
    } else if (n_low < 2L || n_high < 2L) {
      status <- 'insufficient_expression_groups'
    } else if (events < 1L) {
      status <- 'no_events'
    } else {
      logrank_fit <- tryCatch(
        survdiff(Surv(time, event) ~ expression_group, data = survival_data),
        error = function(error) error
      )
      if (inherits(logrank_fit, 'error')) {
        status <- 'logrank_failed'
        warning_text <- conditionMessage(logrank_fit)
      } else {
        logrank_chisq <- unname(logrank_fit$chisq)
        logrank_p <- pchisq(logrank_chisq, df = 1L, lower.tail = FALSE)
      }

      cox_warning <- ''
      cox_fit <- tryCatch(
        withCallingHandlers(
          coxph(Surv(time, event) ~ expression_group, data = survival_data,
                ties = 'efron', x = TRUE),
          warning = function(warning) {
            cox_warning <<- paste(cox_warning, conditionMessage(warning), sep = '; ')
            invokeRestart('muffleWarning')
          }
        ),
        error = function(error) error
      )
      if (!inherits(cox_fit, 'error') && length(coef(cox_fit)) == 1L && is.finite(coef(cox_fit))) {
        hr <- exp(unname(coef(cox_fit)))
        confidence_interval <- suppressMessages(confint(cox_fit))
        hr_lower95 <- exp(confidence_interval[1])
        hr_upper95 <- exp(confidence_interval[2])
        cox_p <- summary(cox_fit)$coefficients[1, 'Pr(>|z|)']
        if (nzchar(cox_warning)) warning_text <- paste(warning_text, cox_warning, sep = '; ')
      } else {
        if (status == 'ok') status <- 'logrank_only'
        cox_message <- if (inherits(cox_fit, 'error')) conditionMessage(cox_fit) else 'nonfinite Cox coefficient'
        warning_text <- paste(warning_text, cox_message, sep = '; ')
      }
    }

    survival_rows[[survival_row_index]] <- data.frame(
      gene_symbol = gene_metadata_de$gene_symbol[gene_index],
      gene_id = gene_metadata_de$gene_id[gene_index],
      stratum = stratum,
      median_expression = median_expression,
      n = nrow(survival_data), n_low = n_low, n_high = n_high, events = events,
      logrank_chisq = logrank_chisq, logrank_p = logrank_p,
      hr = hr, hr_lower95 = hr_lower95, hr_upper95 = hr_upper95,
      cox_p = cox_p, status = status, warning = trimws(warning_text, whitespace = '[; ]'),
      stringsAsFactors = FALSE
    )
  }
  message('Survival completed: ', stratum, ' (', nrow(norm.cts), ' genes)')
}

survival_long <- do.call(rbind, survival_rows)
survival_long$logrank_padj <- NA_real_
survival_long$cox_padj <- NA_real_
for (stratum in survival_strata) {
  stratum_index <- survival_long$stratum == stratum
  finite_logrank <- stratum_index & is.finite(survival_long$logrank_p)
  finite_cox <- stratum_index & is.finite(survival_long$cox_p)
  survival_long$logrank_padj[finite_logrank] <- p.adjust(survival_long$logrank_p[finite_logrank], method = 'BH')
  survival_long$cox_padj[finite_cox] <- p.adjust(survival_long$cox_p[finite_cox], method = 'BH')
}
fwrite(survival_long, file.path(dir_survival, 'survival_all_genes_long.tsv'),
       sep = '\t', na = 'NA')

## 6.3. Convert survival results to one row per gene ####
survival_wide <- data.frame(
  gene_symbol = gene_metadata_de$gene_symbol,
  gene_id = gene_metadata_de$gene_id,
  stringsAsFactors = FALSE
)
for (stratum in survival_strata) {
  stratum_table <- survival_long[survival_long$stratum == stratum, , drop = FALSE]
  stratum_table <- stratum_table[match(survival_wide$gene_id, stratum_table$gene_id), , drop = FALSE]
  stratum_table <- stratum_table[, setdiff(names(stratum_table), c('gene_symbol', 'gene_id', 'stratum')),
                                 drop = FALSE]
  names(stratum_table) <- paste0('survival_', stratum, '_', names(stratum_table))
  survival_wide <- cbind(survival_wide, stratum_table)
}
fwrite(survival_wide, file.path(dir_survival, 'survival_all_genes_wide.tsv'),
       sep = '\t', na = 'NA')

survival_summary <- list(
  endpoint = 'overall survival', event_coding = '0=censored; 1=death', time_unit = 'days',
  expression = 'DESeq2 size-factor-normalized counts',
  cutpoint = 'median within each stratum before survival eligibility filtering',
  split = 'HIGH >= median; LOW < median', hr_contrast = 'HIGH / LOW',
  genes = nrow(norm.cts), samples = ncol(norm.cts), eligible_patients = sum(valid_survival),
  gene_universe = 'exactly the expression-filtered DESeq2 input genes; no DEG significance filter',
  status_counts = as.list(table(survival_long$status))
)
write_json(survival_summary, file.path(dir_survival, 'survival_summary.json'),
           pretty = TRUE, auto_unbox = TRUE, na = 'null')
saveRDS(list(table = survival_long, wide = survival_wide, summary = survival_summary),
        file.path(dir_survival, 'survival.rds'))


# 7. Merge DE and survival results ----------------------------------------
cat("integrated_results\n", file=file.path(output_dir, "progress.log"), append=TRUE)

## 7.1. Confirm the shared gene universe ####
de_input_table <- de_table[de_table$de_filter_pass, , drop = FALSE]
stopifnot(
  !anyDuplicated(de_input_table$gene_id),
  !anyDuplicated(survival_wide$gene_id),
  identical(de_input_table$gene_id, survival_wide$gene_id)
)

## 7.2. Merge by gene ID ####
survival_columns <- setdiff(names(survival_wide), c('gene_symbol', 'gene_id'))
integrated_table <- cbind(
  de_input_table,
  survival_wide[match(de_input_table$gene_id, survival_wide$gene_id), survival_columns, drop = FALSE]
)
fwrite(integrated_table, file.path(dir_merge, 'DE_survival_all_genes.tsv'),
       sep = '\t', na = 'NA')
saveRDS(integrated_table, file.path(dir_merge, 'DE_survival_integrated.rds'))

dim(integrated_table)
head(integrated_table[, intersect(c(
  'gene_symbol', 'gene_id', 'log2FoldChange', 'padj', 'direction',
  'survival_pooled_hr', 'survival_pooled_logrank_padj',
  'survival_TNBC_hr', 'survival_TNBC_logrank_padj'
), names(integrated_table)), drop = FALSE])

## 7.3. Write a portable HTML summary ####
deg_up <- sum(de_table$direction == 'Up_in_TNBC', na.rm = TRUE)
deg_down <- sum(de_table$direction == 'Down_in_TNBC', na.rm = TRUE)
hallmark_significant <- gsea_summary$Hallmark$significant
gobp_significant <- gsea_summary$GOBP$significant

html_report <- c(
  '<!doctype html>', '<html lang="ko"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>TCGA-BRCA TNBC biomarker discovery</title>',
  '<style>body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1100px;margin:40px auto;padding:0 20px;color:#172a40}h1,h2{color:#173f67}.metrics{display:flex;flex-wrap:wrap;gap:12px}.metric{background:#eef4fa;padding:14px 18px;border-radius:8px}img{max-width:100%;border:1px solid #dce3eb;margin:8px 0 24px}a{color:#1964ad}</style>',
  '</head><body>',
  '<h1>TCGA-BRCA · TNBC vs Non-TNBC</h1>',
  '<div class="metrics">',
  sprintf('<div class="metric"><b>Samples</b><br>%s</div>', ncol(cnts)),
  sprintf('<div class="metric"><b>DE/survival genes</b><br>%s</div>', sum(keep)),
  sprintf('<div class="metric"><b>DEGs</b><br>%s up / %s down</div>', deg_up, deg_down),
  sprintf('<div class="metric"><b>Batch model</b><br>%s</div>', selected_batch_id),
  sprintf('<div class="metric"><b>GSEA</b><br>%s Hallmark / %s GOBP</div>', hallmark_significant, gobp_significant),
  '</div>',
  '<h2>QC and batch correction</h2>',
  '<img src="01_qc/before_batch_pca_group.png" alt="PCA before batch correction">',
  '<img src="02_batch/after_batch_pca_group.png" alt="PCA after batch correction">',
  '<p><a href="02_batch/batch_candidates.tsv">Batch candidate table</a></p>',
  '<h2>Differential expression</h2>',
  '<img src="04_de_plots/volcano.png" alt="Volcano plot">',
  if (file.exists(file.path(output_dir, '04_de_plots/DEG_heatmap.png')))
    '<img src="04_de_plots/DEG_heatmap.png" alt="DEG heatmap">' else '<p>DEG heatmap: no significant results to plot.</p>',
  '<p><a href="03_de/DEG.tsv">DEG table</a> · <a href="03_de/DE_all_genes.tsv">All DE results</a></p>',
  '<h2>GSEA</h2>',
  if (file.exists(file.path(output_dir, '05_gsea/Hallmark_NES_top.png')))
    '<img src="05_gsea/Hallmark_NES_top.png" alt="Hallmark GSEA">' else '<p>Hallmark GSEA: no significant results to plot.</p>',
  if (file.exists(file.path(output_dir, '05_gsea/GOBP_NES_top.png')))
    '<img src="05_gsea/GOBP_NES_top.png" alt="GO BP GSEA">' else '<p>GO BP GSEA: no significant results to plot.</p>',
  '<h2>Survival and integrated result</h2>',
  '<p><a href="06_survival/survival_all_genes_long.tsv">Gene-wise survival table</a></p>',
  '<p><a href="07_merge/DE_survival_all_genes.tsv">Integrated DE and survival table</a></p>',
  '<p>Exploratory discovery result. Independent-cohort and single-cell validation are required.</p>',
  '</body></html>'
)
writeLines(html_report, file.path(output_dir, 'index.html'), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(output_dir, 'sessionInfo.txt'))
if (!is.na(script_file)) file.copy(script_file, file.path(output_dir, 'analysis_entrypoint.R'))
report_info <- file.info(file.path(output_dir, 'index.html'))
write_json(list(
  status = 'SUCCEEDED', report = 'index.html', bytes = unname(report_info$size),
  md5 = unname(tools::md5sum(file.path(output_dir, 'index.html'))),
  finished_at = format(Sys.time(), tz = 'UTC', usetz = TRUE)
), file.path(output_dir, 'report_status.json'), pretty = TRUE, auto_unbox = TRUE)

cat('\nAnalysis completed.\n')
cat('Output directory: ', output_dir, '\n', sep = '')
cat('HTML report: ', file.path(output_dir, 'index.html'), '\n', sep = '')
