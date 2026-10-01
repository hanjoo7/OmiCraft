suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1], mustWork=TRUE)
cohort <- fread(file.path(root, 'cohort.tsv'), data.table=FALSE, na.strings=c('', 'NA'))
out <- file.path(root, 'prepared')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_counts <- function(path) fread(path, skip='gene_id', data.table=FALSE)
first <- read_counts(cohort$counts_path[1])
first <- first[grepl('^ENSG', first$gene_id), ]
annotation <- data.frame(gene_id=first$gene_id, gene_symbol=first$gene_name, gene_type=first$gene_type)
counts <- matrix(0L, nrow(first), nrow(cohort), dimnames=list(first$gene_id, cohort$sample_id))
for (i in seq_len(nrow(cohort))) {
  tab <- read_counts(cohort$counts_path[i])
  tab <- tab[match(annotation$gene_id, tab$gene_id), ]
  stopifnot(identical(tab$gene_id, annotation$gene_id), all(is.finite(tab$unstranded)),
            all(tab$unstranded >= 0), all(tab$unstranded == round(tab$unstranded)))
  counts[, i] <- as.integer(tab$unstranded)
  if (i %% 100 == 0) message('Loaded raw counts: ', i, '/', nrow(cohort))
}
metadata <- cohort[, setdiff(names(cohort), c('counts_path','source_md5','file_id')), drop=FALSE]
rownames(metadata) <- metadata$sample_id
clinical <- metadata[, c('patient_id','os_event','os_time_days')]
read_gmt <- function(path, collection) {
  entries <- strsplit(readLines(path), '\t', fixed=TRUE)
  terms <- vapply(entries, function(x) x[1], character(1))
  genes <- lapply(entries, function(x) unique(x[-c(1,2)]))
  list(term2gene=data.frame(term=rep(terms, lengths(genes)), gene=unlist(genes, use.names=FALSE)),
       term2name=data.frame(term=terms, name=terms),
       provenance=list(source='Broad MSigDB', version='2025.1.Hs', collection=collection,
                       file=basename(path), md5=unname(tools::md5sum(path))))
}
gene_sets <- list(Hallmark=read_gmt(file.path(root,'raw/msigdb_hallmark.gmt'),'H'),
                  GOBP=read_gmt(file.path(root,'raw/msigdb_gobp.gmt'),'C5:GO:BP'))
inputs <- list(counts=counts, metadata=metadata, annotation=annotation, clinical=clinical, gene_sets=gene_sets)
for (name in names(inputs)) saveRDS(inputs[[name]], file.path(out, paste0(name,'.rds')), compress=FALSE)
fwrite(annotation, file.path(out, 'annotation.tsv'), sep='\t')
fwrite(metadata, file.path(out, 'metadata.tsv'), sep='\t')
manifest <- list(created_at=format(Sys.time(), tz='UTC', usetz=TRUE), genes=nrow(counts), samples=ncol(counts),
                 groups=as.list(table(metadata$group)), raw_integer_counts=TRUE,
                 files=lapply(names(inputs), function(name) {
                   path <- file.path(out,paste0(name,'.rds'))
                   list(name=name,path=path,md5=unname(tools::md5sum(path)))
                 }))
write_json(manifest,file.path(out,'prepared_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
cat(toJSON(manifest, pretty=TRUE, auto_unbox=TRUE), '\n')
