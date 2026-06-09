###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(bio3d)
  library(optparse)
})
option_list <- list(
  make_option(c('--pdb-dir'), default='structural_complexes'),
  make_option(c('--metadata'), default='metadata/weinstein_numbering_opioids.csv'),
  make_option(c('--out'), default='summary/interhelical_distances.csv')
)
opt <- parse_args(OptionParser(option_list=option_list))

read_bw <- function(path) read.csv(path, sep=';', stringsAsFactors=FALSE, check.names=FALSE)
parse_name <- function(stem) {
  parts <- strsplit(stem, '[-_]')[[1]]
  receptor <- toupper(parts[1])
  partner <- paste(parts[-1], collapse='_')
  cls <- ifelse(grepl('^Arr|^ARR', partner), 'arrestin', ifelse(grepl('^G', partner), 'gprotein', 'other'))
  c(receptor=receptor, partner=partner, partner_class=cls)
}
anchor_res <- function(bw, receptor, tm, bw_decimal) {
  sub <- bw[bw$receptor == receptor & bw$TM == tm,]
  if(nrow(sub) == 0) return(NA_integer_)
  # BW X.50 is stored as x. residue for X.decimal is x + decimal - 50
  as.integer(sub$x[1] + (bw_decimal - 50))
}
ca_xyz <- function(pdb, chain='A', resid) {
  sel <- atom.select(pdb, chain=chain, resno=resid, elety='CA')
  if(length(sel$xyz) < 3) return(c(NA_real_,NA_real_,NA_real_))
  as.numeric(pdb$xyz[1, sel$xyz[1:3]])
}
dist3 <- function(a,b) sqrt(sum((a-b)^2))
angle3 <- function(a,b,c,d) {
  v1 <- b-a; v2 <- d-c
  acos(sum(v1*v2)/(sqrt(sum(v1*v1))*sqrt(sum(v2*v2)))) * 180/pi
}

bw <- read_bw(opt$metadata)
files <- list.files(opt$pdb_dir, pattern='\\.pdb$', full.names=TRUE)
rows <- list()
for(f in files) {
  stem <- tools::file_path_sans_ext(basename(f))
  info <- parse_name(stem)
  rec <- info[['receptor']]
  if(!(rec %in% bw$receptor)) next
  pdb <- tryCatch(read.pdb(f), error=function(e) NULL)
  if(is.null(pdb)) next
  # Anchors: 3.50 as conserved intracellular TM3 anchor; 6.30 and 7.53 capture intracellular opening; H8 optional 8.47.
  r350 <- anchor_res(bw, rec, 'TM3', 50)
  r630 <- anchor_res(bw, rec, 'TM6', 30)
  r753 <- anchor_res(bw, rec, 'TM7', 53)
  r847 <- anchor_res(bw, rec, 'H8', 47)
  # helix vector ends use TM starts/ends from metadata when possible.
  tm3 <- bw[bw$receptor==rec & bw$TM=='TM3',]
  tm6 <- bw[bw$receptor==rec & bw$TM=='TM6',]
  tm7 <- bw[bw$receptor==rec & bw$TM=='TM7',]
  a350 <- ca_xyz(pdb, 'A', r350); a630 <- ca_xyz(pdb, 'A', r630); a753 <- ca_xyz(pdb, 'A', r753); a847 <- ca_xyz(pdb, 'A', r847)
  tm3a <- ca_xyz(pdb, 'A', as.integer(tm3$start[1])); tm3b <- ca_xyz(pdb, 'A', as.integer(tm3$end[1]))
  tm6a <- ca_xyz(pdb, 'A', as.integer(tm6$start[1])); tm6b <- ca_xyz(pdb, 'A', as.integer(tm6$end[1]))
  tm7a <- ca_xyz(pdb, 'A', as.integer(tm7$start[1])); tm7b <- ca_xyz(pdb, 'A', as.integer(tm7$end[1]))
  rows[[length(rows)+1]] <- data.frame(
    complex=stem, receptor=rec, partner=info[['partner']], partner_class=info[['partner_class']],
    tm3_anchor=r350, tm6_anchor=r630, tm7_anchor=r753, h8_anchor=r847,
    tm3_tm6=dist3(a350,a630), tm3_tm7=dist3(a350,a753), tm7_h8=dist3(a753,a847),
    tm3_tm6_angle=angle3(tm3a,tm3b,tm6a,tm6b), tm3_tm7_angle=angle3(tm3a,tm3b,tm7a,tm7b),
    stringsAsFactors=FALSE
  )
}
out <- if(length(rows)) do.call(rbind, rows) else data.frame()
dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
write.csv(out, opt$out, row.names=FALSE)
cat('Wrote', opt$out, 'with', nrow(out), 'rows\n')
