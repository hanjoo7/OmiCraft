suppressPackageStartupMessages({library(Matrix);library(data.table)})
root <- commandArgs(trailingOnly=TRUE)[1]
expr <- as(readMM(file.path(root,'expression.mtx')), 'CsparseMatrix')
meta <- fread(file.path(root,'cells.tsv'),data.table=FALSE)
genes <- fread(file.path(root,'genes.tsv'),data.table=FALSE)
stopifnot(nrow(expr)==nrow(genes),ncol(expr)==nrow(meta),!anyDuplicated(genes$feature_name),
          all(is.finite(expr@x)),all(expr@x>=0),all(expr@x==round(expr@x)))
rownames(expr) <- genes$feature_name
colnames(expr) <- meta$cell_id
rownames(meta) <- meta$cell_id
saveRDS(expr,file.path(root,'expression.rds'),compress=FALSE)
saveRDS(meta,file.path(root,'metadata.rds'),compress=FALSE)
cat(nrow(expr),'genes x',ncol(expr),'cells converted\n')
