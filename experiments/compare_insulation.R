#!/usr/bin/env Rscript
# Compare Crane insulation at Level-1 TAD boundaries: SpectralTAD vs SpectralTAD_TV.
# Same input, window, and k-selection; only the affinity (contacts vs TV-ADMM C) differs.

suppressPackageStartupMessages({
  library(SpectralTAD)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list(
    sparse = "/nfs/turbo/umms-minjilab/lpullela/559/559_spectralTAD/data/lateG1_chr22_10kb_NONE.sparse.txt",
    chr = "chr22",
    resolution = 10000L,
    window_size = 200L,
    delta = 25L,
    n_perm = 1000L,
    tv_max_iter = 50L,
    tv_gamma = 0.1,
    out_dir = "/nfs/turbo/umms-minjilab/lpullela/559/559_spectralTAD/results"
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    out[[gsub("-", "_", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  out$resolution <- as.integer(out$resolution)
  out$window_size <- as.integer(out$window_size)
  out$delta <- as.integer(out$delta)
  out$n_perm <- as.integer(out$n_perm)
  out$tv_max_iter <- as.integer(out$tv_max_iter)
  out$tv_gamma <- as.numeric(out$tv_gamma)
  out
}

opt <- parse_args(args)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(0)

message("Reading ", opt$sparse)
raw <- as.matrix(read.table(opt$sparse, header = FALSE, sep = "\t",
                            colClasses = c("numeric", "numeric", "numeric")))
colnames(raw) <- c("bin1", "bin2", "counts")
M <- HiCcompare::sparse2full(raw)
storage.mode(M) <- "double"
n <- nrow(M)
coords <- as.numeric(colnames(M))
message("Matrix: ", n, " x ", n)

call_tads <- function(fun, label, extra = list()) {
  message("\n=== ", label, " ===")
  t0 <- proc.time()[[3]]
  args <- c(list(cont_mat = raw, chr = opt$chr, levels = 1L, qual_filter = FALSE,
                 z_clust = FALSE, resolution = opt$resolution,
                 window_size = opt$window_size), extra)
  res <- do.call(fun, args)
  elapsed <- proc.time()[[3]] - t0
  bed <- as.data.frame(res$Level_1)
  message(label, ": ", nrow(bed), " TADs in ", round(elapsed, 1), "s")
  attr(bed, "elapsed") <- elapsed
  bed
}

bed_ref <- call_tads(SpectralTAD, "SpectralTAD")
bed_tv  <- call_tads(SpectralTAD_TV, "SpectralTAD_TV",
                     list(tv_gamma = opt$tv_gamma, tv_max_iter = opt$tv_max_iter))

write.table(bed_ref, file.path(opt$out_dir, "insulation_SpectralTAD_Level_1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bed_tv, file.path(opt$out_dir, "insulation_SpectralTAD_TV_Level_1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Crane 2015 insulation: contacts across bin i in a delta x delta square.
insulation <- function(M, delta) {
  n <- nrow(M)
  ins <- rep(NA_real_, n)
  for (i in seq.int(delta + 1L, n - delta)) {
    ins[i] <- sum(M[(i - delta):(i - 1L), i:(i + delta - 1L)])
  }
  mu <- mean(ins, na.rm = TRUE)
  log2(pmax(ins, 1e-12) / mu)
}

# Internal TAD ends mapped onto matrix bin indices (boundary sits at the start of the next TAD).
boundary_idx <- function(bed, coords, resolution) {
  b <- sort(unique(bed$end[-nrow(bed)]))
  idx <- match(b, coords)
  idx[is.finite(idx)]
}

is_local_min <- function(ins, idx) {
  ok <- idx > 1L & idx < length(ins)
  idx <- idx[ok]
  ins[idx] < ins[idx - 1L] & ins[idx] < ins[idx + 1L]
}

eval_bounds <- function(ins, idx, delta, n_perm) {
  n <- length(ins)
  lo <- delta + 1L
  hi <- n - delta
  idx <- idx[idx >= lo & idx <= hi]
  obs <- mean(ins[idx], na.rm = TRUE)
  frac <- mean(is_local_min(ins, idx), na.rm = TRUE)
  L <- hi - lo + 1L
  rel <- idx - lo
  null_means <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    s <- sample.int(L, 1L) - 1L
    nb <- ((rel + s) %% L) + lo
    null_means[p] <- mean(ins[nb], na.rm = TRUE)
  }
  z <- (obs - mean(null_means)) / sd(null_means)
  list(n = length(idx), mean = obs, z = z, frac = frac,
       null_mean = mean(null_means), null_sd = sd(null_means))
}

message("\nComputing insulation (delta = ", opt$delta, " bins)")
ins <- insulation(M, opt$delta)

idx_ref <- boundary_idx(bed_ref, coords, opt$resolution)
idx_tv  <- boundary_idx(bed_tv, coords, opt$resolution)
m_ref <- eval_bounds(ins, idx_ref, opt$delta, opt$n_perm)
m_tv  <- eval_bounds(ins, idx_tv,  opt$delta, opt$n_perm)

fmt <- function(m) {
  sprintf("n=%d  mean=%.4f  z=%.2f  frac_local_min=%.3f", m$n, m$mean, m$z, m$frac)
}

lines <- c(
  paste0("lateG1 chr22  ", opt$resolution, " bp  window=", opt$window_size,
         "  insulation_delta=", opt$delta, "  n_perm=", opt$n_perm),
  paste0("SpectralTAD     ", fmt(m_ref), "  time_s=", round(attr(bed_ref, "elapsed"), 1)),
  paste0("SpectralTAD_TV  ", fmt(m_tv),  "  time_s=", round(attr(bed_tv, "elapsed"), 1),
         "  gamma=", opt$tv_gamma, "  max_iter=", opt$tv_max_iter),
  "More negative z is better (deeper insulation than a spacing-preserving circular-shift null).",
  "frac_local_min is the share of boundaries sitting at an insulation local minimum."
)
writeLines(lines, file.path(opt$out_dir, "insulation_compare_chr22.txt"))
cat(paste(lines, collapse = "\n"), "\n")
