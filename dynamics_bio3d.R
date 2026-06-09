###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

#!/usr/bin/env Rscript
# Bio3D backend for the reconstructed opioid receptor dynamics pipeline.
#
# This is the reference backend because the original 2021 pipeline performed
# Normal Mode Analysis in R/Bio3D. It recomputes domain-level fluctuation
# fold-changes and Bhattacharyya coefficients from the PDB files, without using
# the legacy .RData objects. The .RData files are only retained for validation
# and for reproducing the published legacy plots.

suppressPackageStartupMessages({
  library(bio3d)
  library(optparse)
})

option_list <- list(
  make_option(c("--root"), type="character", default=".", help="Repository root"),
  make_option(c("--cutoff"), type="double", default=10.0, help="Reserved for parity with Python backend"),
  make_option(c("--limit"), type="integer", default=NA, help="Debug: process first N complexes"),
  make_option(c("--outdir"), type="character", default="processed_results/dynamics_bio3d", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
root <- normalizePath(opt$root, mustWork=TRUE)
outdir <- file.path(root, opt$outdir)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(root, "summary"), recursive=TRUE, showWarnings=FALSE)

domain_order <- c("TM1","ICL1","TM2","ECL1","TM3","ICL2","TM4","ECL2","TM5","ICL3","TM6","ECL3","TM7","H8","Other")

read_meta <- function(file) {
  read.csv(file.path(root, "metadata", file), sep=";", stringsAsFactors=FALSE, check.names=FALSE)
}

bw <- read_meta("weinstein_numbering_opioids.csv")
bw$start <- as.integer(bw$start); bw$end <- as.integer(bw$end); bw$x <- as.integer(bw$x)

parse_complex <- function(stem) {
  parts <- strsplit(stem, "[-_]")[[1]]
  receptor <- toupper(parts[1])
  partner <- sub(paste0("^", receptor, "[-_]"), "", stem, ignore.case=TRUE)
  class <- ifelse(grepl("^ARR", partner, ignore.case=TRUE), "arrestin", ifelse(grepl("^G", partner, ignore.case=TRUE), "gprotein", "other"))
  list(receptor=receptor, partner=partner, partner_class=class)
}

infer_template <- function(partner) {
  p <- toupper(partner)
  if (grepl("6PWC", p)) return("6PWC")
  if (grepl("6U1N", p)) return("6U1N")
  return("6DDF")
}

monomer_for <- function(receptor, partner) {
  dyn <- file.path(root, "dynamics_complexes")
  candidates <- c(
    file.path(dyn, paste0(receptor, "_", infer_template(partner), ".pdb")),
    file.path(dyn, paste0(receptor, "_6DDF.pdb")),
    file.path(dyn, paste0(receptor, "_6PWC.pdb")),
    file.path(dyn, paste0(receptor, "_6U1N.pdb"))
  )
  candidates[file.exists(candidates)][1]
}

domain_for <- function(receptor, resno) {
  sub <- bw[toupper(bw$receptor) == toupper(receptor), ]
  for (i in seq_len(nrow(sub))) {
    if (!is.na(sub$start[i]) && resno >= sub$start[i] && resno <= sub$end[i]) return(as.character(sub$TM[i]))
  }
  sub <- sub[order(sub$start),]
  if (nrow(sub) >= 2) {
    for (i in seq_len(nrow(sub)-1)) {
      a <- as.character(sub$TM[i]); b <- as.character(sub$TM[i+1])
      if (resno > sub$end[i] && resno < sub$start[i+1]) {
        if (a == "TM1" && b == "TM2") return("ICL1")
        if (a == "TM3" && b == "TM4") return("ICL2")
        if (a == "TM5" && b == "TM6") return("ICL3")
        return(paste0(a,"-",b," loop"))
      }
    }
    if (resno > sub$end[nrow(sub)]) return("H8")
  }
  "Other"
}

safe_fluct <- function(pdb_file) {
  pdb <- read.pdb(pdb_file)
  ca <- atom.select(pdb, chain="A", elety="CA")
  if (length(ca$atom) < 10) stop("fewer than 10 CA atoms in chain A")
  pdb_ca <- trim.pdb(pdb, inds=ca)
  # Bio3D NMA. nma.pdb is the historical Bio3D route used for coarse-grained NMA.
  modes <- nma.pdb(pdb_ca, ff="calpha")
  fl <- fluct.nma(modes)
  atoms <- pdb_ca$atom
  data.frame(resno=as.integer(atoms$resno), resid=atoms$resid, fluct=as.numeric(fl), stringsAsFactors=FALSE)
}

normalize <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) == 0) return(numeric())
  x <- x - min(x, na.rm=TRUE)
  if (sum(x, na.rm=TRUE) <= 0) x <- rep(1, length(x))
  x / sum(x, na.rm=TRUE)
}

bc_coef <- function(p, q) {
  p <- normalize(p); q <- normalize(q); n <- min(length(p), length(q))
  if (n == 0) return(NA_real_)
  sum(sqrt(p[seq_len(n)] * q[seq_len(n)]))
}

complex_files <- list.files(file.path(root, "dynamics_complexes"), pattern=".*-.*\\.pdb$", full.names=TRUE)
complex_files <- sort(complex_files)
if (!is.na(opt$limit)) complex_files <- head(complex_files, opt$limit)

residue_rows <- list(); domain_rows <- list(); i <- 1; j <- 1
for (cf in complex_files) {
  stem <- sub("\\.pdb$", "", basename(cf))
  info <- parse_complex(stem)
  if (!(info$receptor %in% c("DOR","KOR","MOR","NOP"))) next
  mf <- monomer_for(info$receptor, info$partner)
  if (length(mf) == 0 || is.na(mf)) { cat("No monomer for", stem, "\n"); next }
  cat("Bio3D dynamics:", stem, "vs", basename(mf), "\n")
  ok <- tryCatch({
    comp <- safe_fluct(cf); mono <- safe_fluct(mf)
    common <- intersect(comp$resno, mono$resno)
    comp <- comp[match(common, comp$resno), ]; mono <- mono[match(common, mono$resno), ]
    domains <- vapply(common, function(r) domain_for(info$receptor, r), character(1))
    fold <- comp$fluct / mono$fluct
    residue_rows[[i]] <- data.frame(
      complex=stem, monomer=sub("\\.pdb$", "", basename(mf)), receptor=info$receptor,
      partner=info$partner, partner_class=info$partner_class, resno=common,
      residue=comp$resid, domain=domains, fluctuation_complex=comp$fluct,
      fluctuation_monomer=mono$fluct, fold_change=fold, stringsAsFactors=FALSE
    ); i <- i + 1
    for (dom in unique(domains)) {
      idx <- which(domains == dom)
      domain_rows[[j]] <- data.frame(
        complex=stem, receptor=info$receptor, partner=info$partner, partner_class=info$partner_class,
        domain=dom, n_residues=length(idx), mean_fluctuation_complex=mean(comp$fluct[idx], na.rm=TRUE),
        mean_fluctuation_monomer=mean(mono$fluct[idx], na.rm=TRUE), mean_fold_change=mean(fold[idx], na.rm=TRUE),
        bhattacharyya_coefficient=bc_coef(comp$fluct[idx], mono$fluct[idx]), stringsAsFactors=FALSE
      ); j <- j + 1
    }
    TRUE
  }, error=function(e) { cat("WARNING Bio3D failed for", stem, ":", conditionMessage(e), "\n"); FALSE })
}

residue_df <- if (length(residue_rows)) do.call(rbind, residue_rows) else data.frame()
domain_df <- if (length(domain_rows)) do.call(rbind, domain_rows) else data.frame()
write.csv(residue_df, file.path(outdir, "bio3d_residue_fluctuations.csv"), row.names=FALSE)
write.csv(domain_df, file.path(outdir, "bio3d_domain_dynamics_metrics.csv"), row.names=FALSE)

# Convert Bio3D results to the same legacy-style wide tables.
if (nrow(domain_df) > 0) {
  domains_keep <- c("TM1","TM2","TM3","TM4","TM5","TM6","TM7","H8")
  for (cls in c("gprotein", "arrestin")) {
    sub <- domain_df[domain_df$partner_class == cls & domain_df$domain %in% domains_keep, ]
    if (nrow(sub) == 0) next
    complexes <- unique(sub$complex)
    fluc <- data.frame(Complex=complexes, stringsAsFactors=FALSE)
    bc <- data.frame(row.names=complexes)
    for (dom in domains_keep) {
      vals_f <- tapply(sub$mean_fold_change[sub$domain == dom], sub$complex[sub$domain == dom], mean)
      vals_b <- tapply(sub$bhattacharyya_coefficient[sub$domain == dom], sub$complex[sub$domain == dom], mean)
      fluc[[paste0("Fluc_", dom)]] <- as.numeric(vals_f[complexes])
      bc[[dom]] <- as.numeric(vals_b[complexes])
    }
    fluc$DXR <- sub("[-_].*$", "", fluc$Complex)
    fluc$Partner <- sub("^[A-Za-z0-9]+[-_]", "", fluc$Complex)
    prefix <- ifelse(cls == "gprotein", "gprot", "arrestin")
    write.csv(fluc, file.path(outdir, paste0(prefix, "_fluc_fold_change_bio3d.csv")), row.names=FALSE)
    write.csv(bc, file.path(outdir, paste0(prefix, "_bc_scores_bio3d.csv")))
  }
}
cat("Bio3D dynamics complete. Outputs in", outdir, "\n")
