#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/inspect_rdata.R file.RData [outdir]")
rdata <- args[[1]]
outdir <- ifelse(length(args) >= 2, args[[2]], "processed_results/rdata_inspection")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(rdata, envir=e)
objs <- ls(e)
base <- tools::file_path_sans_ext(basename(rdata))
summary <- data.frame(file=basename(rdata), object=character(), class=character(), typeof=character(), length=integer(), nrow=integer(), ncol=integer(), stringsAsFactors=FALSE)
for (nm in objs) {
  obj <- get(nm, envir=e)
  cl <- paste(class(obj), collapse="|")
  nr <- ifelse(is.null(dim(obj)), NA, dim(obj)[1])
  nc <- ifelse(is.null(dim(obj)) || length(dim(obj)) < 2, NA, dim(obj)[2])
  summary <- rbind(summary, data.frame(file=basename(rdata), object=nm, class=cl, typeof=typeof(obj), length=length(obj), nrow=nr, ncol=nc, stringsAsFactors=FALSE))
  sink(file.path(outdir, paste0(base, "__", nm, "__str.txt")))
  cat("FILE:", rdata, "\nOBJECT:", nm, "\nCLASS:", cl, "\n\n")
  print(utils::str(obj, max.level=3))
  sink()
}
write.csv(summary, file=file.path(outdir, paste0(base, "__object_summary.csv")), row.names=FALSE)
cat("Inspected", rdata, "->", outdir, "\n")
