#!/usr/bin/env Rscript
# Run SpectralTAD on a 3-column sparse contact matrix (bin start, bin start, count).

suppressPackageStartupMessages({
  library(SpectralTAD)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list(
    sparse = NULL,
    chr = NULL,
    resolution = 10000,
    levels = 2,
    window_size = 25,
    out_dir = "results"
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    val <- args[[i + 1]]
    out[[gsub("-", "_", key)]] <- val
    i <- i + 2
  }
  out$resolution <- as.integer(out$resolution)
  out$levels <- as.integer(out$levels)
  out$window_size <- as.integer(out$window_size)
  out
}

opt <- parse_args(args)
if (is.null(opt$sparse) || is.null(opt$chr)) {
  stop("Usage: run_spectraltad.R --sparse FILE --chr CHR [--resolution 10000] [--levels 2] [--window-size 25] [--out-dir DIR]")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading sparse matrix: ", opt$sparse)
mat <- as.matrix(read.table(opt$sparse, header = FALSE, sep = "\t", colClasses = c("numeric", "numeric", "numeric")))
colnames(mat) <- c("bin1", "bin2", "counts")
message("Contacts: ", nrow(mat))

# At 10 kb, 25 bins is only 250 kb; 200 bins ~ 2 Mb, the default TAD-size window.
if (opt$resolution <= 10000 && opt$window_size == 25) {
  opt$window_size <- as.integer(2000000 / opt$resolution)
  message("Adjusted window_size to ", opt$window_size, " bins (~2 Mb)")
}

chr_label <- if (startsWith(opt$chr, "chr")) opt$chr else paste0("chr", opt$chr)
out_prefix <- file.path(opt$out_dir, paste0(chr_label, "_", opt$resolution))

message("Calling SpectralTAD on ", chr_label, " at ", opt$resolution, " bp")
tads <- SpectralTAD(
  cont_mat = mat,
  chr = chr_label,
  levels = opt$levels,
  qual_filter = FALSE,
  z_clust = FALSE,
  resolution = opt$resolution,
  window_size = opt$window_size,
  out_format = "bedpe",
  out_path = paste0(out_prefix, "_tads")
)

for (nm in names(tads)) {
  f <- paste0(out_prefix, "_", nm, ".tsv")
  write.table(tads[[nm]], f, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", f, " (", nrow(tads[[nm]]), " TADs)")
}

print(lapply(tads, head))
