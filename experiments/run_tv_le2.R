#!/usr/bin/env Rscript
# SpectralTAD_TV with γ=0.1, λz=0.1, λe=2.0, μ=1, σ=1; insulation vs existing SpectralTAD.

suppressPackageStartupMessages(library(SpectralTAD))

sparse  <- "/nfs/turbo/umms-minjilab/lpullela/559/559_spectralTAD/data/lateG1_chr22_10kb_NONE.sparse.txt"
out_dir <- "/nfs/turbo/umms-minjilab/lpullela/559/559_spectralTAD/results"
delta   <- 25L
n_perm  <- 1000L
set.seed(0)

raw <- as.matrix(read.table(sparse, header = FALSE, sep = "\t",
                            colClasses = c("numeric", "numeric", "numeric")))
colnames(raw) <- c("bin1", "bin2", "counts")
M <- HiCcompare::sparse2full(raw)
storage.mode(M) <- "double"
n <- nrow(M)
coords <- as.numeric(colnames(M))
message("Matrix: ", n, " x ", n)

message("=== SpectralTAD_TV  gamma=0.1  lambda_z=0.1  lambda_e=2  mu=1  sigma=1  (λ / σ1(Y)) ===")
t0 <- proc.time()[[3]]
tv <- SpectralTAD_TV(
  cont_mat = raw, chr = "chr22", levels = 1L, qual_filter = FALSE, z_clust = FALSE,
  resolution = 10000L, window_size = 200L,
  tv_gamma = 0.1, tv_lambda_z = 0.1, tv_lambda_e = 2.0, tv_mu = 1.0, tv_sigma = 1.0,
  tv_max_iter = 50L
)
elapsed <- proc.time()[[3]] - t0
bed_tv <- as.data.frame(tv$Level_1)
message("SpectralTAD_TV: ", nrow(bed_tv), " TADs in ", round(elapsed, 1), "s")
write.table(bed_tv, file.path(out_dir, "insulation_SpectralTAD_TV_sigma1_Level_1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bed_ref <- read.table(file.path(out_dir, "insulation_SpectralTAD_Level_1.tsv"),
                      header = TRUE, sep = "\t")
bed_unscaled <- read.table(file.path(out_dir, "insulation_SpectralTAD_TV_le2_Level_1.tsv"),
                           header = TRUE, sep = "\t")

raw_ins <- rep(NA_real_, n)
for (i in seq.int(delta + 1L, n - delta)) {
  raw_ins[i] <- sum(M[(i - delta):(i - 1L), i:(i + delta - 1L)])
}
covered <- which(is.finite(raw_ins) & raw_ins > 0)
ins <- rep(NA_real_, n)
ins[covered] <- log2(raw_ins[covered] / mean(raw_ins[covered]))

boundary_idx <- function(bed) {
  b <- sort(unique(bed$end[-nrow(bed)]))
  idx <- match(b, coords)
  idx[is.finite(idx)]
}
is_local_min <- function(idx) {
  idx <- idx[idx > 1L & idx < n]
  ins[idx] < ins[idx - 1L] & ins[idx] < ins[idx + 1L]
}
eval_bounds <- function(idx) {
  idx <- intersect(idx, covered)
  obs <- mean(ins[idx])
  frac <- mean(is_local_min(idx), na.rm = TRUE)
  pos <- match(idx, covered)
  L <- length(covered)
  null_means <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    s <- sample.int(L, 1L) - 1L
    nb <- covered[((pos + s - 1L) %% L) + 1L]
    null_means[p] <- mean(ins[nb])
  }
  z <- (obs - mean(null_means)) / sd(null_means)
  list(n = length(idx), mean = obs, z = z, frac = frac, null_mean = mean(null_means))
}

m_ref <- eval_bounds(boundary_idx(bed_ref))
m_old <- eval_bounds(boundary_idx(bed_unscaled))
m_new <- eval_bounds(boundary_idx(bed_tv))

fmt <- function(lab, m, n_tads) {
  sprintf("%-36s TADs=%-3d  bounds=%-3d  ins_mean=%.4f  z=%.2f  frac_local_min=%.3f",
          lab, n_tads, m$n, m$mean, m$z, m$frac)
}

lines <- c(
  "lateG1 chr22  10000 bp  window=200  insulation_delta=25  n_perm=1000",
  "TV params: γ=0.1  λz=0.1  λe=2  μ=1  σ=1  max_iter=50",
  "Null: circular shift along coverage-masked bins. More negative z is better.",
  fmt("SpectralTAD", m_ref, nrow(bed_ref)),
  fmt("SpectralTAD_TV  unscaled λ", m_old, nrow(bed_unscaled)),
  fmt("SpectralTAD_TV  λ / σ1(Y)", m_new, nrow(bed_tv)),
  sprintf("TV σ1-scaled runtime: %.1fs", elapsed)
)
writeLines(lines, file.path(out_dir, "insulation_compare_chr22_sigma1.txt"))
cat(paste(lines, collapse = "\n"), "\n")
