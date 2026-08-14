#' SpectralTAD_TV: SpectralTAD with TV-ADMM affinity
#'
#' Identical to SpectralTAD() except that the normalized-Laplacian affinity
#' is built from the TV-ADMM coefficient matrix C (from ssc_admm_nuc_tv in
#' the companion Python script) rather than from the raw contact submatrix.
#'
#' After solving for C, the symmetric affinity W = |C| + |C^T| is formed and
#' normalized as D^{-1/2} W D^{-1/2} (same Laplacian step as the original).
#' All downstream logic — eigenvector computation, unit-circle projection,
#' gap distances, silhouette/z-score cut selection — is unchanged.
#'
#' @param tv_gamma     TV regularisation weight gamma (default 0.1).
#' @param tv_lambda_e  Weight on ||E||_1 before 1/σ₁(Y) scaling (default 1.0).
#' @param tv_lambda_z  Weight on reconstruction loss before 1/σ₁(Y) scaling (default 0.1).
#' @param tv_mu        ADMM penalty for X = C_off constraint (default 1.0).
#' @param tv_sigma     ADMM penalty for TV auxiliary constraints (default 1.0).
#' @param tv_max_iter  Maximum ADMM iterations (default 50).
#' @param tv_tol       Convergence tolerance on Frobenius residuals (default 1e-4).
#'
#' All other parameters are identical to SpectralTAD().
#'
#' @export
#' @import dplyr
#' @import magrittr
#' @importFrom GenomicRanges GRanges GRangesList start
#' @importFrom utils write.table

# ── TV-ADMM helpers (R port of ssc_admm_nuc_total_var.py) ────────────────────

.soft_thresh_tv = function(x, tau) sign(x) * pmax(abs(x) - tau, 0)

.finite_diff_mat = function(N) {
  D = matrix(0, N - 1, N)
  idx = seq_len(N - 1)
  D[cbind(idx, idx)]     = -1
  D[cbind(idx, idx + 1)] =  1
  D
}

# Returns the N x N coefficient matrix C from ssc_admm_nuc_tv.
# Y is the (n x N) contact window (n = N for a symmetric Hi-C submatrix).
.tv_admm_C = function(Y,
                      lambda_e = 1.0,
                      lambda_z = 0.1,
                      gamma    = 0.1,
                      mu       = 1.0,
                      sigma    = 1.0,
                      max_iter = 50L,
                      tol      = 1e-4) {
  N = ncol(Y)
  D = .finite_diff_mat(N)
  K = crossprod(D)                                # t(D) %*% D,  N x N
  eig_K = eigen(K, symmetric = TRUE)
  V     = eig_K$vectors
  ev    = eig_K$values                            # ascending from eigen()
  denom = mu + sigma * outer(ev, ev, "+")         # N x N
  A_inv = solve(lambda_z * crossprod(Y) + mu * diag(N))

  X = matrix(0, N, N); C = matrix(0, N, N); E = matrix(0, nrow(Y), N)
  P = matrix(0, N - 1, N); Q = matrix(0, N, N - 1)
  Lambda = matrix(0, N, N); Pi_P = matrix(0, N - 1, N); Pi_Q = matrix(0, N, N - 1)

  for (it in seq_len(max_iter)) {
    X_prev = X
    C_off  = C; diag(C_off) = 0

    # 1. X-update
    X = A_inv %*% (lambda_z * (t(Y) %*% (Y - E)) + mu * C_off - Lambda)

    # 2. C-update (Sylvester equation solved via eigenbasis of K)
    P_tilde = P - Pi_P / sigma
    Q_tilde = Q - Pi_Q / sigma
    RHS_C   = mu * (X + Lambda / mu) + sigma * (t(D) %*% P_tilde + Q_tilde %*% D)
    C = V %*% ((t(V) %*% RHS_C %*% V) / denom) %*% t(V)
    diag(C) = 0

    # 3-4. P- and Q-updates
    DC  = D %*% C;  CDt = C %*% t(D)
    P   = .soft_thresh_tv(DC  + Pi_P / sigma, gamma / sigma)
    Q   = .soft_thresh_tv(CDt + Pi_Q / sigma, gamma / sigma)

    # 5. E-update
    E = .soft_thresh_tv(Y - Y %*% X, lambda_e / lambda_z)

    # 6. Dual updates
    C_off  = C; diag(C_off) = 0
    Lambda = Lambda + mu    * (X   - C_off)
    Pi_P   = Pi_P   + sigma * (DC  - P)
    Pi_Q   = Pi_Q   + sigma * (CDt - Q)

    # Convergence
    res = max(norm(X - C_off, "F"), norm(DC - P, "F"), norm(CDt - Q, "F"))
    if (res < tol && mu * norm(X - X_prev, "F") < tol) break
  }
  C
}

# ── Main function (mirrors SpectralTAD) ──────────────────────────────────────

SpectralTAD_TV = function(cont_mat, chr, levels = 1, qual_filter = FALSE,
                          z_clust = FALSE, eigenvalues = 2, min_size = 5,
                          window_size = 25,
                          resolution = "auto", gap_threshold = 1,
                          grange = FALSE, out_format = "none", out_path = chr,
                          tv_gamma    = 0.1,
                          tv_lambda_e = 1.0,
                          tv_lambda_z = 0.1,
                          tv_mu       = 1.0,
                          tv_sigma    = 1.0,
                          tv_max_iter = 50L,
                          tv_tol      = 1e-4) {

  options(scipen = 999)

  if (missing("chr")) {
    stop("Must specify chromosome")
  }

  row_test = dim(cont_mat)[1]
  col_test = dim(cont_mat)[2]

  if (row_test == col_test) {
    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }
  }

  if (col_test == 3) {

    if (!is.matrix(cont_mat)) {
      cont_mat = as.matrix(cont_mat)
    }

    message("Converting to n x n matrix")

    if (nrow(cont_mat) == 1) {
      stop("Matrix is too small to convert to full")
    }
    cont_mat = HiCcompare::sparse2full(cont_mat)

    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }

    if (resolution == "auto") {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat))-dplyr::lag(as.numeric(colnames(cont_mat)))))[1])
    }

  } else if (col_test-row_test == 3) {

    message("Converting to n x n matrix")

    start_coords = cont_mat[,2]
    resolution   = as.numeric(cont_mat[1,3])-as.numeric(cont_mat[1,2])
    cont_mat     = as.matrix(cont_mat[,-c(seq_len(3))])

    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }

    colnames(cont_mat) = start_coords

  } else if (col_test!=3 & (row_test != col_test) & (col_test-row_test != 3)) {

    stop("Contact matrix must be sparse or n x n or n x (n+3)!")

  } else if ( (resolution == "auto") & (col_test-row_test == 0) ) {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat))-dplyr::lag(as.numeric(colnames(cont_mat)))))[1])
  }

  if (resolution>200000) {
    stop("Resolution must be less than (or equal to) 200kb")
  }

  if (nrow(cont_mat) < 2000000/resolution) {
    stop("Matrix must be larger than 2 megabases divided by resolution")
  }

  bed = .windowedSpec_TV(cont_mat, chr = chr, resolution = resolution,
                         z_clust = z_clust, eigenvalues = eigenvalues,
                         min_size = min_size, window_size = window_size,
                         qual_filter = qual_filter, gap_threshold = gap_threshold,
                         tv_gamma = tv_gamma, tv_lambda_e = tv_lambda_e,
                         tv_lambda_z = tv_lambda_z, tv_mu = tv_mu,
                         tv_sigma = tv_sigma, tv_max_iter = tv_max_iter,
                         tv_tol = tv_tol) %>%
    mutate(Level = 1)

  coords = cbind(match(bed$start, as.numeric(colnames(cont_mat))), match(bed$end-resolution, as.numeric(colnames(cont_mat))))

  tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])

  called_tads = list(bed)

  curr_lev = 2

  while (curr_lev != (levels + 1) ) {

    coords = cbind(match(called_tads[[curr_lev-1]]$start, as.numeric(colnames(cont_mat))), match(called_tads[[curr_lev-1]]$end-resolution, as.numeric(colnames(cont_mat))))

    less_5 = which( (coords[,2]-coords[,1])<min_size*2  )

    if (length(less_5)>0) {
      pres_tads = called_tads[[curr_lev-1]][less_5,]
      coords    = coords[-less_5, ]
    } else {
      pres_tads = c()
    }

    if (is.null(nrow(coords))) {
      coords = t(as.matrix(coords))
    }

    tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])

    zeros = which(unlist(lapply(tads, function(x) nrow(x)-sum(rowSums(x)==0)))<min_size*2)

    if (length(zeros)>0) {
      pres_tads = rbind(pres_tads, called_tads[[curr_lev-1]][zeros,])
      tads[zeros] = NULL
    }

    sub_tads = lapply(tads, function(x) {
      .windowedSpec_TV(x, chr = chr, resolution = resolution, qual_filter = qual_filter,
                       z_clust = TRUE, min_size = min_size,
                       tv_gamma = tv_gamma, tv_lambda_e = tv_lambda_e,
                       tv_lambda_z = tv_lambda_z, tv_mu = tv_mu,
                       tv_sigma = tv_sigma, tv_max_iter = tv_max_iter,
                       tv_tol = tv_tol)
    })

    called_tads[[curr_lev]] = bind_rows(sub_tads, pres_tads) %>% mutate(Level = curr_lev) %>% arrange(start)

    curr_lev = curr_lev+1
  }

  names(called_tads) = paste0("Level_", seq_len(levels))

  if ( !(out_format == "none")) {
    if (out_format %in% c("bedpe", "juicebox")) {
      bed_out = bind_rows(called_tads) %>%
        dplyr::select(chr,start,end)
      bed_out1 <- bed_out
      colnames(bed_out1) <- c("chr1", "start1", "end1")
      bed_out = bind_cols(bed_out, bed_out1) %>%
        mutate(name = ".", score = ".", strand1 =".", strand2 = ".")
      bound_tads = bind_rows(called_tads)
      colors = c("0,0,0", "255,0,0", "0,255,0", "0,0,255")
      bed_out = bed_out %>%
        mutate(color =colors[bound_tads$Level])
      write.table(bed_out, out_path, quote = FALSE,
                  row.names = FALSE, sep = "\t", col.names = FALSE)
    } else if (out_format %in% c("bed", "hicexplorer")) {
      bed_out = bind_rows(called_tads) %>%
        dplyr::select(chr,start,end)
      write.table(bed_out, out_path, quote = FALSE,
                  row.names = FALSE, sep = "\t", col.names = FALSE)
    } else {
      warning("No file output, unsupported output format chosen")
    }
  }

  if (grange == TRUE) {
    called_tads = lapply(called_tads, function(x) {
      GenomicRanges::GRanges(x)
    })
    called_tads = GenomicRanges::GRangesList(called_tads)
  }

  return(called_tads)
}


# ── Internal sliding-window function ─────────────────────────────────────────

.windowedSpec_TV = function(cont_mat, resolution, chr,
                            gap_filter = TRUE, z_clust = FALSE, qual_filter = TRUE,
                            eigenvalues = 2, min_size = 5,
                            window_size = ceiling(2000000/resolution),
                            gap_threshold = 1,
                            tv_gamma    = 0.1,
                            tv_lambda_e = 1.0,
                            tv_lambda_z = 0.1,
                            tv_mu       = 1.0,
                            tv_sigma    = 1.0,
                            tv_max_iter = 50L,
                            tv_tol      = 1e-4) {

  window_size = ceiling(window_size)

  Group_over = dplyr::bind_rows()

  start    = 1
  end      = window_size
  end_loop = 0

  if (end+window_size>nrow(cont_mat)) {
    end = nrow(cont_mat)
  }

  while (end_loop == 0) {

    sub_filt = cont_mat[seq(start,end, 1), seq(start,end, 1)]

    zero_thresh   = round(nrow(sub_filt)*(gap_threshold))
    non_gaps_within = which((colSums(sub_filt == 0))<zero_thresh)

    sub_filt = sub_filt[non_gaps_within, non_gaps_within]

    if (length(nrow(sub_filt)) == 0) {
      start = end
      end   = start+window_size
      if (start == nrow(cont_mat)) { end_loop = 1; next }
      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) { end = nrow(cont_mat) }
      next
    }

    if (nrow(sub_filt) < min_size*2) {
      start = end
      end   = start+window_size
      if (start == nrow(cont_mat)) { end_loop = 1; next }
      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) { end = nrow(cont_mat) }
      next
    }

    dist_sub = 1/(1+sub_filt)

    # ── TV-ADMM substitution ────────────────────────────────────────────────
    # Run the TV-ADMM solver on the raw contact window. λ_e and λ_z are divided
    # by σ₁(Y) so the E soft-threshold is scale-invariant across windows
    # (same as tad_common.tv_affinity). Then symmetrise C and build the
    # normalised Laplacian from it instead of from sub_filt.

    sigma1 = svd(sub_filt, nu = 0, nv = 0)$d[1]
    if (!is.finite(sigma1) || sigma1 < 1e-10) {
      start = end
      end = start + window_size
      if (start == nrow(cont_mat)) { end_loop = 1; next }
      if ((end + (2000000 / resolution)) > nrow(cont_mat)) { end = nrow(cont_mat) }
      next
    }

    C_tv = .tv_admm_C(sub_filt,
                      lambda_e = tv_lambda_e / sigma1,
                      lambda_z = tv_lambda_z / sigma1,
                      gamma    = tv_gamma,
                      mu       = tv_mu,
                      sigma    = tv_sigma,
                      max_iter = tv_max_iter,
                      tol      = tv_tol)

    W  = abs(C_tv) + abs(t(C_tv))          # symmetric affinity

    dr = rowSums(W)
    Dinvsqrt = diag(1 / sqrt(pmax(dr, 1e-12)))   # guard zero-degree rows
    sub_mat  = Dinvsqrt %*% W %*% Dinvsqrt

    # ── End substitution; rest of .windowedSpec is unchanged ────────────────

    colnames(sub_mat) = colnames(cont_mat)[non_gaps_within]
    sub_mat[is.nan(sub_mat)] = 0

    Eigen = get_eigs(sub_mat, NEig = eigenvalues)

    eig_vals = Eigen$values
    eig_vecs = Eigen$vectors

    large_small = order(-eig_vals)
    eig_vals    = eig_vals[large_small]
    eig_vecs    = eig_vecs[,large_small]

    index     = 1
    Group_mem = list()

    clusters = seq_len(ceiling( (end-start+1)/min_size))

    norm_ones = sqrt(dim(sub_mat)[2])

    for (i in seq_len(dim(eig_vecs)[2])) {
      eig_vecs[,i] = (eig_vecs[,i]/sqrt(sum(eig_vecs[,i]^2)))  * norm_ones
      if (eig_vecs[1,i] !=0) {
        eig_vecs[,i] = -1*eig_vecs[,i] * sign(eig_vecs[1,i])
      }
    }

    n = dim(eig_vecs)[1]
    k = dim(eig_vecs)[2]

    eig_vecs = crossprod(diag(diag(tcrossprod(eig_vecs))^(-1/2)), eig_vecs)

    point_dist = sqrt(rowSums( (eig_vecs-rbind(NA,eig_vecs[-nrow(eig_vecs),]))^2  ))

    if (z_clust) {

      sig_bounds = which(scale(point_dist[-length(point_dist)])>2)
      sig_bounds = subset(sig_bounds, sig_bounds>min_size)
      dist_bounds = which(c(min_size*2,diff(sig_bounds))<min_size)

      if (length(dist_bounds) > 0) {
        sig_bounds = sig_bounds[-dist_bounds]
      }

      TAD_start = c(1, sig_bounds+1)
      TAD_end   = c(sig_bounds, nrow(sub_filt))
      widths    = (TAD_end-TAD_start)+1
      memberships = unlist(lapply(seq_len(length(TAD_start)), function(x) rep(x,widths[x])))

      if (length(sig_bounds) == 0) {
        end_group = dplyr::bind_rows()
      } else {
        sig_bounds  = which(scale(point_dist[-length(point_dist)])>2)
        sig_bounds  = subset(sig_bounds, sig_bounds>min_size)
        dist_bounds = which(c(min_size*2,diff(sig_bounds))<min_size)
        end_group   = data.frame(ID = as.numeric(colnames(sub_filt)), Group = memberships)
        end_group   = end_group %>% dplyr::mutate(group_place = Group) %>% dplyr::group_by(group_place) %>% dplyr::mutate(Group = last(ID)) %>% dplyr::ungroup() %>% dplyr::select(ID, Group)
      }

    } else {

      gap_order = order(-point_dist)
      sil_score = c()

      for (cluster in clusters) {

        k = 1
        partition_found = 0
        first_run       = TRUE
        cutpoints       = c()

        while(partition_found == 0) {

          new_gap   = gap_order[k]
          cutpoints = c(cutpoints, new_gap)

          diff_points = which( abs(new_gap-cutpoints[-length(cutpoints)]) <= min_size)

          if (length(diff_points)>0) {
            cutpoints = cutpoints[-length(cutpoints)]
          }

          if (length(cutpoints) == cluster) {
            partition_found = 1
          } else {
            k = k+1
          }
        }

        if (any(is.na(cutpoints))) { next }

        cutpoints  = cutpoints[order(cutpoints)]
        cutpoints  = c(1, cutpoints, length(non_gaps_within)+1)
        group_size = diff(cutpoints)

        memberships = c()
        for (i in seq_len(length(group_size))) {
          memberships = c(memberships, rep(i,times = group_size[i]))
        }

        sil = summary(cluster::silhouette(memberships,dist_sub))
        sil_score = c(sil_score, sil$si.summary[4])
        Group_mem[[cluster]] = memberships
      }

      end_group = Group_mem[[which(diff(sil_score)<0)[1]]]

      if (length(end_group) == 0) {
        end_group = dplyr::bind_rows()
      } else {
        end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = end_group)
        end_group = end_group %>%dplyr::mutate(group_place = Group) %>%dplyr::group_by(group_place) %>%dplyr::mutate(Group = max(ID)) %>% ungroup() %>% dplyr::select(ID, Group)
      }
    }

    if (end == nrow(cont_mat)) {
      Group_over = dplyr::bind_rows(Group_over, end_group)
      end_loop   = 1
    } else {

      if (nrow(end_group)!=0) {
        end_IDs = which(end_group$Group == last(end_group$Group))
      } else {
        end_IDs = 1:window_size
      }

      start = end-length(end_IDs)+1

      if (length(start) == 0 ) {
        start = end
      }

      if (nrow(end_group != 0)) {
        end = start+window_size
      } else {
        end = start+window_size*2
      }

      end_group  = end_group[-end_IDs, ]
      Group_over = dplyr::bind_rows(Group_over, end_group)

      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) {
        end = nrow(cont_mat)
      }
    }
  }

  if (z_clust) {
    if (nrow(Group_over) > 0) {
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>%dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
        dplyr::filter((end-start)/resolution >= min_size) %>%dplyr::arrange(start)
    } else {
      bed = Group_over
    }
  } else {
    if (qual_filter) {
      fin_range    = match(Group_over$ID,colnames(cont_mat))
      over_dist_mat = 1/(1+cont_mat[fin_range, fin_range])
      sil          = cluster::silhouette(Group_over$Group, over_dist_mat)
      ave_sil      = summary(sil)$clus.avg.widths
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
        dplyr::mutate(Sil_Score = ave_sil) %>% dplyr::filter( ((end-start)/resolution >= min_size) & Sil_Score > .15)  %>%dplyr::arrange(start)
    } else {
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%dplyr::filter((end-start)/resolution >= min_size) %>% dplyr::arrange(start)
    }
  }

  return(bed)
}
