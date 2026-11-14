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
arrays <- read.celfiles(files, verbose=FALSE)
matrix <- exprs(rma(arrays, target="core"))
stopifnot(identical(colnames(matrix), basename(files)), all(is.finite(matrix)))
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
