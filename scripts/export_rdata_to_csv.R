#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript scripts/export_rdata_to_csv.R file.RData outdir")
rdata <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(rdata, envir=e)
base <- tools::file_path_sans_ext(basename(rdata))
domains <- c("TM1","TM2","TM3","TM4","TM5","TM6","TM7","H8")

write_obj <- function(obj, prefix, outdir){
  if (is.matrix(obj) || is.array(obj)) {
    df <- as.data.frame(obj, stringsAsFactors=FALSE)
    if (ncol(df) == 8 && grepl("bc_scores", prefix)) colnames(df) <- paste0("BC_", domains)
    if (!is.null(rownames(obj))) df <- cbind(Complex=rownames(obj), df)
    fn <- file.path(outdir, paste0(prefix, ".csv"))
    write.csv(df, fn, row.names=FALSE)
    return(fn)
  }
  if (is.data.frame(obj)) {
    df <- obj
    if (!is.null(rownames(obj)) && !all(rownames(obj) == as.character(seq_len(nrow(obj))))) {
      df <- cbind(RowName=rownames(obj), df)
    }
    fn <- file.path(outdir, paste0(prefix, ".csv"))
    write.csv(df, fn, row.names=FALSE)
    return(fn)
  }
  if (is.vector(obj) && !is.list(obj)) {
    fn <- file.path(outdir, paste0(prefix, ".csv"))
    write.csv(data.frame(index=seq_along(obj), value=obj), fn, row.names=FALSE)
    return(fn)
  }
  return(character())
}

manifest <- data.frame(source=character(), object=character(), class=character(), exported_file=character(), stringsAsFactors=FALSE)
for (nm in ls(e)) {
  obj <- get(nm, envir=e)
  safe <- gsub("[^A-Za-z0-9_.-]", "_", nm)
  files <- write_obj(obj, paste0(base, "__", safe), outdir)
  if (length(files)>0) {
    manifest <- rbind(manifest, data.frame(source=basename(rdata), object=nm, class=paste(class(obj), collapse="|"), exported_file=files, stringsAsFactors=FALSE))
  }
}
write.csv(manifest, file=file.path(outdir, paste0(base, "__export_manifest.csv")), row.names=FALSE)
cat("Exported", nrow(manifest), "objects from", rdata, "to", outdir, "\n")
