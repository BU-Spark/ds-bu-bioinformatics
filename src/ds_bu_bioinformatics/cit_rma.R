# Normalize the full deposited assay batch; preserve transcript-cluster IDs.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(pd.hugene.2.0.st))
suppressPackageStartupMessages(library(hugene20sttranscriptcluster.db))
pins <- c(oligo="1.74.0", pd.hugene.2.0.st="3.14.1",
          hugene20sttranscriptcluster.db="8.8.0", AnnotationDbi="1.72.0",
          org.Hs.eg.db="3.22.0", affxparser="1.82.0", preprocessCore="1.72.0")
observed <- vapply(names(pins), function(p) as.character(packageVersion(p)), character(1))
stopifnot(getRversion()=="4.5.3", identical(observed, pins))
files <- readLines(args[1L])
stopifnot(length(files)>1L, !anyDuplicated(basename(files)))
qc <- data.frame(cel_file=basename(files), cel_qc_pass=TRUE, cel_qc_reason="readable_finite_nonnegative_raw_signals")
arrays <- tryCatch(read.celfiles(files, verbose=FALSE), error=function(e) NULL)
if (is.null(arrays)) {
    gc()
    # Diagnose individual source files without repairing, imputing or converting
    # damaged raw measurements. Keep all original files in the source cache.
    for (i in seq_along(files)) {
        problem <- tryCatch({
            one <- read.celfiles(files[i], verbose=FALSE)
            stopifnot(all(is.finite(exprs(one))), all(exprs(one)>=0))
            NULL
        }, error=function(e) conditionMessage(e))
        if (!is.null(problem)) {
            qc$cel_qc_pass[i] <- FALSE
            qc$cel_qc_reason[i] <- paste("unreadable_or_invalid_CEL:",problem)
            message("Excluded raw array ", basename(files[i]), ": ", problem)
        }
        if (exists("one",inherits=FALSE)) rm(one)
        if (i %% 10L == 0L) gc()
    }
    gc()
    stopifnot(sum(qc$cel_qc_pass)>1L)
    arrays <- read.celfiles(files[qc$cel_qc_pass], verbose=FALSE)
}
bad <- !vapply(seq_len(ncol(arrays)),function(i) all(is.finite(exprs(arrays)[,i])) && all(exprs(arrays)[,i]>=0),logical(1))
if (any(bad)) {
    ids <- match(sampleNames(arrays)[bad],qc$cel_file)
    stopifnot(!anyNA(ids))
    qc$cel_qc_pass[ids] <- FALSE
    qc$cel_qc_reason[ids] <- "nonfinite_or_negative_raw_CEL_signal"
    arrays <- arrays[,!bad]
}
write.table(qc,file.path(args[2L],"cel_qc.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
stopifnot(ncol(arrays)>1L,identical(sampleNames(arrays),basename(files)[qc$cel_qc_pass]))
matrix <- exprs(rma(arrays, target="core"))
stopifnot(identical(colnames(matrix), basename(files)[qc$cel_qc_pass]), all(is.finite(matrix)))
con <- file(file.path(args[2L],"expression.tsv"),"w")
writeLines(paste(c("feature_id",colnames(matrix)),collapse="\t"),con)
for (i in seq_len(nrow(matrix))) {
    writeLines(paste(c(rownames(matrix)[i],sprintf("%.17g",matrix[i,])),collapse="\t"),con)
}
close(con)
annotations <- AnnotationDbi::select(hugene20sttranscriptcluster.db, rownames(matrix),
                                     c("ENTREZID","SYMBOL"), "PROBEID")
write.table(annotations,file.path(args[2L],"annotations.tsv"),sep="\t",quote=FALSE,
            row.names=FALSE,na="")
writeLines(c(R.version.string, paste(names(observed), observed)),file.path(args[2L],"versions.txt"))
