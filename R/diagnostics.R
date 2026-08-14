## -----------------------------------------------------------------------------
#' Initialize incremental target-clade counts
#'
#' @keywords internal
#' @noRd
initialize_target_tree_asdsf_state <- function(phy_target,
                                                edge_list_chains,
                                                rooted = TRUE) {
    counts <- lapply(edge_list_chains, function(chain_edges) {
        if (!length(chain_edges)) return(rep(0, phy_target$Nnode))
        prop_clades_par(
            E_target = phy_target$edge,
            edges = chain_edges,
            rooted = rooted,
            normalize = FALSE
        )
    })
    list(
        counts = do.call(rbind, counts),
        totals = lengths(edge_list_chains)
    )
}

#' Update incremental target-clade counts
#' @keywords internal
#' @noRd
update_target_tree_asdsf_state <- function(state, phy_target,
                                            new_edge_list_chains,
                                            rooted = TRUE) {
    for (chain_id in seq_along(new_edge_list_chains)) {
        chain_edges <- new_edge_list_chains[[chain_id]]
        if (!length(chain_edges)) next
        state$counts[chain_id, ] <- state$counts[chain_id, ] + prop_clades_par(
            E_target = phy_target$edge,
            edges = chain_edges,
            rooted = rooted,
            normalize = FALSE
        )
        state$totals[[chain_id]] <- state$totals[[chain_id]] + length(chain_edges)
    }
    state
}

#' Compute ASDSF from incremental target-clade counts
#' @keywords internal
#' @noRd
target_tree_asdsf_from_state <- function(state, min_freq = 0) {
    freq_matrix <- state$counts / pmax(state$totals, 1L)

    if (nrow(freq_matrix) <= 1L) {
        per_clade_sd <- rep(0, ncol(freq_matrix))
    } else {
        per_clade_sd <- apply(freq_matrix, 2, stats::sd)
        per_clade_sd[is.na(per_clade_sd)] <- 0
    }

    keep_mask <- apply(freq_matrix, 2, function(vals) any(vals > min_freq))
    if (!any(keep_mask)) 0 else mean(per_clade_sd[keep_mask])
}

#' Compute ASDSF across chains using clades from a target tree
#'
#' @param phy_target A rooted `phylo` object defining the reference clades.
#' @param edge_list_chains List of chains, each a list of edge matrices.
#' @param rooted Logical; treat trees as rooted when matching clades.
#' @param normalize Logical; compute normalized clade frequencies. Only `TRUE`
#'   is supported.
#' @param ncores Number of threads used for clade counting.
#' @param min_freq Minimum clade frequency used when averaging SDs.
#' @return The average standard deviation of target-clade frequencies across
#'   chains.
#' @keywords internal
compute_target_tree_asdsf <- function(phy_target,
                                      edge_list_chains,
                                      rooted = TRUE,
                                      normalize = TRUE,
									  ncores = 1,
                                      min_freq = 0) {

    if (!inherits(phy_target, "phylo")) {
        stop('`phy_target` must be a phylo object')
    }
    if (!is.list(edge_list_chains) || length(edge_list_chains) == 0) {
        stop('`edge_list_chains` must be a non-empty list of chains')
    }

    if (!isTRUE(normalize)) stop('ASDSF requires normalized clade frequencies')

	RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    state <- initialize_target_tree_asdsf_state(
        phy_target, edge_list_chains, rooted = rooted)
    target_tree_asdsf_from_state(state, min_freq = min_freq)
}

#' Compute clade retention curve from tree
#' Returns data frame with confidence values and cumulative count of clades remaining
#' @keywords internal
#' @noRd
compute_retention_curve <- function(tree) {
    # Get confidence values from node labels
    conf_values <- as.numeric(tree$node.label)
    conf_values <- conf_values[!is.na(conf_values)]
    total_clades <- length(conf_values)

    # Sort descending to get unique thresholds
    conf_sorted <- sort(unique(conf_values), decreasing = TRUE)

    # For each unique confidence value, count how many clades remain at that threshold
    n_remaining <- vapply(conf_sorted, function(thresh) {
        sum(conf_values >= thresh)
    }, numeric(1))

    data.frame(
        confidence = conf_sorted,
        n_remaining = n_remaining,
        frac_remaining = n_remaining / total_clades
    )
}

#' Plot retention curves
#'
#' @param retention_data A single retention data frame, or a named list of them.
#'   When a list, names are used as legend labels.
#' @param sample_name Plot title.
#' @param cutoff Optional vertical line at a confidence threshold.
#' @keywords internal
#' @noRd
plot_retention_curve <- function(retention_data, sample_name = '', cutoff = NULL) {

    if (is.data.frame(retention_data)) {
        retention_data <- list(retention_data)
    }
    if (is.null(names(retention_data)) || any(names(retention_data) == "")) {
        base_name <- if (!is.null(sample_name) && nzchar(sample_name)) sample_name else "retention"
        if (length(retention_data) == 1L) {
            names(retention_data) <- base_name
        } else {
            names(retention_data) <- paste0(base_name, "_", seq_along(retention_data))
        }
    }
    df <- do.call(rbind, lapply(names(retention_data), function(nm) {
        d <- retention_data[[nm]]
        if (!all(c("confidence", "n_remaining") %in% colnames(d))) {
            stop("retention_data entries must contain columns: confidence, n_remaining")
        }
        d$group <- nm
        d
    }))

    single <- length(retention_data) == 1L

    if (single) {
        p <- ggplot(df, aes(x = confidence, y = n_remaining)) +
            geom_step(direction = 'hv', linewidth = 0.8, color = 'steelblue')
    } else {
        p <- ggplot(df, aes(x = confidence, y = n_remaining, color = group)) +
            geom_step(direction = 'hv', linewidth = 0.8)
    }

    if (!is.null(cutoff)) {
        cutoff <- as.numeric(cutoff)
        cutoff <- cutoff[!is.na(cutoff)]
        if (length(cutoff) == 0L) {
            stop("cutoff must contain at least one numeric value")
        }
        p <- p + geom_vline(xintercept = cutoff, color = 'red', linetype = 'dashed', linewidth = 0.6)
    }

    p <- p +
        scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
        labs(title = sample_name, x = 'Confidence threshold', y = 'Clades remaining', color = NULL) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            legend.position = if (single) 'none' else 'bottom',
            panel.grid.minor = element_blank()
        )

    return(p)
}


#' Build a clade matrix from a phylo object
#' @param tree A phylo object
#' @param tip_order A character vector of tip labels
#' @param include_root A logical indicating whether to include the root clade
#' @param add_node_names If TRUE, add and use node labels for row names (default TRUE);
#'   if FALSE, use numeric clade IDs.
#' @return A sparse matrix of clade x tip incidence
#' @keywords internal
#' @noRd
build_clade_matrix <- function(tree, tip_order, include_root = FALSE, add_node_names = TRUE) {
	stopifnot(inherits(tree, "phylo"))
	
	N <- length(tree$tip.label)

	# map original tip indices -> column positions in tip_order
	pos <- rep.int(NA_integer_, N)
	idx_keep <- match(tip_order, tree$tip.label)
	pos[idx_keep] <- seq_along(idx_keep)

	# internal nodes' tip sets (original indices)
	int_nodes <- (N + 1):(N + tree$Nnode)
	int_sets <- lapply(int_nodes, function(v) phangorn::Descendants(tree, v, type = "tips")[[1]])
	int_sets2 <- lapply(int_sets, function(S) { x <- pos[S]; x[!is.na(x)] })

	if (!include_root) {
		root_mask <- lengths(int_sets2) == length(tip_order)
        int_sets2 <- int_sets2[!root_mask]
        int_nodes <- int_nodes[!root_mask]
	}

	# enforce minor clades: for each set, if size > half of all tips, take its complement
	n_all <- length(tip_order)
	all_idx <- seq_len(n_all)
	sets_minor <- lapply(int_sets2, function(S) {
		if (length(S) > n_all / 2) setdiff(all_idx, S) else S
	})

	# Use filtered internal nodes
	sets <- sets_minor
  if (add_node_names) {
    tree <- add_node_names(tree, prefix = "Node")
    row_names <- tree$node.label[int_nodes - N]
  } else {
    row_names <- as.character(int_nodes)
  }

	# build sparse node x tip incidence
	i <- rep.int(seq_along(sets), lengths(sets))
	j <- unlist(sets, use.names = FALSE)
	M <- Matrix::sparseMatrix(
		i = i, j = j, x = 1L,
		dims = c(length(sets), length(tip_order)),
		dimnames = list(row_names, tip_order)
	)
	
	return(M)
}



#' Build a clade matrix from a cell x variant VAF matrix
#'
#' Treat each variant as defining a clade consisting of cells whose VAF is
#' greater than or equal to a minimum cutoff.
#'
#' @param vaf_mat A numeric matrix (or data.frame coercible to matrix) of
#'   cells x variants. Row names must be cell IDs and column names must be
#'   variant IDs.
#' @param tip_order Optional character vector of cell IDs to define column
#'   order. If NULL, uses rownames(vaf_mat).
#' @param min_vaf Minimum VAF cutoff to include a cell in a variant clade.
#'   Defaults to 0.
#' @param min_cells Minimum number of cells required to be included in a
#'   variant clade (i.e., VAF > min_vaf) for the variant to be retained.
#'   Defaults to 2.
#'
#' @return A sparse matrix of clade (variant) x tip (cell) incidence.
#' @keywords internal
#' @noRd
build_variant_clade_matrix <- function(vaf_mat, tip_order = NULL, min_vaf = 0, min_cells = 2) {

  if (is.data.frame(vaf_mat)) vaf_mat <- as.matrix(vaf_mat)
  stopifnot(is.matrix(vaf_mat))
  if (is.null(rownames(vaf_mat))) stop("vaf_mat must have rownames (cell IDs).")
  if (is.null(colnames(vaf_mat))) stop("vaf_mat must have colnames (variant IDs).")
  if (!is.numeric(min_vaf) || length(min_vaf) != 1L || is.na(min_vaf)) {
    stop("min_vaf must be a single non-NA numeric value.")
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1L || is.na(min_cells) || min_cells < 1) {
    stop("min_cells must be a single numeric value >= 1.")
  }
 
  if (is.null(tip_order)) {
    tip_order <- rownames(vaf_mat)
  } else {
    missing_tips <- setdiff(tip_order, rownames(vaf_mat))
    if (length(missing_tips) > 0L) {
      stop("tip_order contains cells not present in vaf_mat: ", paste(missing_tips, collapse = ", "))
    }
  }

  X <- vaf_mat[tip_order, , drop = FALSE]
  storage.mode(X) <- "numeric"
  X[is.na(X)] <- -Inf

  # pass: cells x variants logical
  pass <- X > min_vaf

  # Keep only variants with at least min_cells passing
  keep_vars <- colSums(pass) >= min_cells
  if (!any(keep_vars)) {
    M <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = integer(0),
      dims = c(0, length(tip_order)),
      dimnames = list(character(0), tip_order)
    )
    return(M)
  }
  pass <- pass[, keep_vars, drop = FALSE]
  kept_variants <- colnames(X)[keep_vars]

  idx <- which(pass, arr.ind = TRUE)

  # Build sparse variant x cell incidence
  M <- Matrix::sparseMatrix(
    i = if (nrow(idx) == 0L) integer(0) else idx[, 2],
    j = if (nrow(idx) == 0L) integer(0) else idx[, 1],
    x = 1L,
    dims = c(length(kept_variants), length(tip_order)),
    dimnames = list(kept_variants, tip_order)
  )

  return(M)
}


#' Compute a Jaccard matrix between variant-defined clades and tree clades
#'
#' Builds a clade x tip incidence matrix from a variants x cells VAF matrix
#' (one clade per variant, cells included if VAF > min_vaf), and compares it to
#' the clades induced by a phylogenetic tree on the shared set of cells.
#'
#' @param vmat A numeric matrix (or data.frame coercible to matrix) of
#'   variants x cells. Row names must be variant IDs and column names must be
#'   cell IDs.
#' @param phy A phylo object.
#' @param min_vaf Minimum VAF cutoff used to include a cell in a variant clade.
#'   Defaults to 0.
#'
#' @return A matrix of Jaccard scores with rows = variants and cols = tree clades.
#' @keywords internal
#' @noRd
variant_tree_jaccard_matrix <- function(vmat, phy, min_vaf = 0, min_cells = 2) {
  stopifnot(inherits(phy, "phylo"))
  if (is.data.frame(vmat)) vmat <- as.matrix(vmat)
  stopifnot(is.matrix(vmat))

  # shared tips (stable order)
  common <- sort(intersect(colnames(vmat), phy$tip.label))
  if (length(common) == 0L) stop("vmat and tree share no cell IDs.")

  phy <- phy %>% ape::keep.tip(common)

  # Variant clades: variants x cells incidence based on VAF > min_vaf,
  # with variants filtered to have >= 2 cells passing min_vaf.
  clade_mat_var <- build_variant_clade_matrix(t(vmat), tip_order = common, min_vaf = min_vaf, min_cells = min_cells)

  clade_mat_tree <- build_clade_matrix(phy, tip_order = common)

  if (nrow(clade_mat_var) == 0L || nrow(clade_mat_tree) == 0L) {
    out <- matrix(0, nrow = nrow(clade_mat_var), ncol = nrow(clade_mat_tree))
    rownames(out) <- rownames(clade_mat_var)
    colnames(out) <- rownames(clade_mat_tree)
    return(out)
  }

  # clade sizes
  var_sizes <- as.numeric(Matrix::rowSums(clade_mat_var))
  tree_sizes <- as.numeric(Matrix::rowSums(clade_mat_tree))

  # intersections: variant x tree-clade counts (sparse)
  inter <- clade_mat_var %*% Matrix::t(clade_mat_tree)  # dgCMatrix
  intersection_sizes <- as.matrix(inter)

  # Union = size1 + size2 - intersection
  union_sizes <- outer(var_sizes, tree_sizes, `+`) - intersection_sizes

  # Jaccard = intersection / union
  J <- intersection_sizes / union_sizes
  J[union_sizes == 0] <- 0

  rownames(J) <- rownames(clade_mat_var)
  colnames(J) <- rownames(clade_mat_tree)
  return(J)
}


#' Compute variant-to-tree precision/recall/F1 curve across confidence cutoffs
#'
#' Construct a variant VAF matrix from long-format mutation counts, compare
#' variant partitions against tree clades using Jaccard similarity, and sweep
#' node-confidence cutoffs to summarize precision, recall, and F1.
#'
#' @param phy_annot A rooted \code{phylo} object with clade confidence stored in
#'   \code{node.label}.
#' @param mut_dat Long-format mutation table containing at least columns
#'   \code{cell}, \code{variant}, \code{a}, and \code{d}.
#' @param min_vaf Minimum VAF used when defining variant-positive cells
#'   (default \code{0.01}).
#' @param min_cells Minimum number of variant-positive cells required for a
#'   variant to be evaluated (default \code{2}).
#' @param j_thres Jaccard threshold for declaring a match between a variant
#'   split and a tree clade (default \code{0.66}).
#' @param n_points Number of confidence cutoffs to evaluate (default
#'   \code{50}).
#' @param ncores Number of cores for sweeping confidence cutoffs. Uses
#'   \code{parallel::mclapply} when \code{ncores > 1} on non-Windows systems.
#'
#' @return A data frame with columns \code{conf_cutoff}, \code{n_pred},
#'   \code{n_tp}, \code{n_tru}, \code{n_recover}, \code{precision},
#'   \code{recall}, \code{f1}, and \code{fraction_well_matched}.
#'
#' @export
compute_variant_pr_curve <- function(
    phy_annot,
    mut_dat,
    min_vaf = 0.01,
    min_cells = 2,
    j_thres = 0.5,
    n_points = 50,
    ncores = 1
) {
    if (!inherits(phy_annot, "phylo")) stop("phy_annot must be a 'phylo' object")

    mut_dat <- mut_dat
    mut_dat$vaf <- mut_dat$a / mut_dat$d
    mut_dat$vaf[is.na(mut_dat$vaf) | is.infinite(mut_dat$vaf)] <- NA_real_

    # Provided by mitodrift: tools/mitodrift/R/utils.R
    vmat <- long_to_mat(mut_dat, "vaf")

    common_cells <- intersect(colnames(vmat), phy_annot$tip.label)
    common_cells <- sort(common_cells)
    if (length(common_cells) == 0L) stop("No overlapping cells between mut_dat and phy_annot tip labels.")

    vmat <- vmat[, common_cells, drop = FALSE]
    tr <- ape::keep.tip(phy_annot, common_cells)

    node_conf <- suppressWarnings(as.numeric(tr$node.label))
    node_conf[is.na(node_conf)] <- 0
    node_ids <- (length(tr$tip.label) + 1):(length(tr$tip.label) + tr$Nnode)
    names(node_conf) <- as.character(node_ids)

    J_full <- variant_tree_jaccard_matrix(
        vmat,
        tr,
        min_vaf = min_vaf,
        min_cells = min_cells
    )

    conf_all <- stats::setNames(node_conf, paste0("Node", names(node_conf)))

    df <- compute_pr_from_jaccard_mat(
        J_full,
        conf_all,
        j_thres = j_thres,
        n_points = n_points,
        ncores = ncores
    )

    var_max_j <- apply(J_full, 1, max)
    df$fraction_well_matched <- mean(var_max_j >= j_thres)
    df
}

#' Compute precision/recall/F1 curve from a truth tree and an estimated tree
#'
#' @param phy_tru A rooted \code{phylo} object representing the ground-truth tree.
#' @param phy_est A rooted \code{phylo} object with clade confidence in \code{node.label}.
#' @param j_thres Jaccard threshold for declaring a match (default 0.5).
#' @param n_points Number of confidence cutoff points to evaluate (default 50).
#' @param ncores Number of cores for sweeping confidence cutoffs. Uses
#'   \code{parallel::mclapply} when \code{ncores > 1} on non-Windows systems.
#'
#' @return Data frame with columns: conf_cutoff, n_pred, n_tp, n_tru,
#'   n_recover, precision, recall, f1.
#' @keywords internal
#' @noRd
compute_tree_pr_curve <- function(phy_tru, phy_est, j_thres = 0.5, n_points = 50, ncores = 1) {
    if (!inherits(phy_tru, "phylo")) stop("phy_tru must be a 'phylo' object")
    if (!inherits(phy_est, "phylo")) stop("phy_est must be a 'phylo' object")

    common_tips <- sort(intersect(phy_tru$tip.label, phy_est$tip.label))
    if (length(common_tips) == 0L) stop("Trees share no tip labels.")

    phy_tru <- ape::keep.tip(phy_tru, common_tips)
    phy_est <- ape::keep.tip(phy_est, common_tips)

    J <- clade_jaccard_matrix(phy_tru, phy_est, add_node_names = FALSE)

    # Map column node IDs to confidence values
    n_tips <- length(phy_est$tip.label)
    col_ids <- as.integer(colnames(J))
    conf_all <- stats::setNames(
        suppressWarnings(as.numeric(phy_est$node.label[col_ids - n_tips])),
        colnames(J)
    )
    conf_all[is.na(conf_all)] <- 0

    compute_pr_from_jaccard_mat(J, conf_all, j_thres = j_thres, n_points = n_points, ncores = ncores)
}

#' Compute precision/recall/F1 curve from a Jaccard matrix and clade confidences
#'
#' @param J Matrix of Jaccard scores (rows = truth clades, cols = predicted clades).
#' @param conf Named numeric vector of clade confidence values keyed by node ID.
#'   Names must match column names of \code{J}.
#' @param j_thres Jaccard threshold for declaring a match (default 0.66).
#' @param n_points Number of confidence cutoff points to evaluate (default 50).
#' @param ncores Number of cores for sweeping confidence cutoffs. Uses
#'   \code{parallel::mclapply} when \code{ncores > 1} on non-Windows systems.
#'
#' @return Data frame with columns: conf_cutoff, n_pred, n_tp, n_tru,
#'   n_recover, precision, recall, f1.
#' @keywords internal
#' @noRd
compute_pr_from_jaccard_mat <- function(J, conf, j_thres = 0.5, n_points = 50, ncores = 1) {

    conf_cols <- conf[colnames(J)]
    conf_cols[is.na(conf_cols)] <- 0

    # Per predicted clade: TP if best-matching truth clade >= j_thres
    col_supported <- apply(J, 2, max) >= j_thres

    # Per truth clade: max confidence among predicted clades that match it
    row_support_conf <- apply(J, 1, function(x) {
        idx <- which(x >= j_thres)
        if (length(idx) == 0L) return(-Inf)
        max(conf_cols[idx], na.rm = TRUE)
    })

    # Sweep confidence cutoffs
    probs <- seq(0, 1, length.out = n_points)
    conf_cutoffs <- as.numeric(stats::quantile(conf_cols, probs = probs, na.rm = TRUE))
    conf_cutoffs[1] <- 0
    conf_cutoffs[length(conf_cutoffs)] <- 1
    conf_cutoffs <- pmin(pmax(conf_cutoffs, 0), 1)

    n_tru <- nrow(J)
    score_for_cutoff <- function(ct) {
        keep <- conf_cols >= ct
        n_pred <- sum(keep)
        n_tp <- sum(col_supported & keep)
        n_recover <- sum(row_support_conf >= ct)
        precision <- if (n_pred == 0L) NA_real_ else n_tp / n_pred
        recall <- if (n_tru == 0L) NA_real_ else n_recover / n_tru
        c(n_pred = n_pred, n_tp = n_tp, n_tru = n_tru,
          n_recover = n_recover, precision = precision, recall = recall)
    }

    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) ncores <- 1L

    score_list <- parallel::mclapply(conf_cutoffs, score_for_cutoff, mc.cores = ncores)
    score_mat <- do.call(cbind, score_list)

    df <- data.frame(
        conf_cutoff = conf_cutoffs,
        n_pred = score_mat["n_pred", ],
        n_tp = score_mat["n_tp", ],
        n_tru = score_mat["n_tru", ],
        n_recover = score_mat["n_recover", ],
        precision = score_mat["precision", ],
        recall = score_mat["recall", ]
    )
    denom <- df$precision + df$recall
    df$f1 <- ifelse(is.na(denom) | denom == 0, NA_real_,
                    2 * df$precision * df$recall / denom)
    df
}

#' Plot precision/recall/F1 versus confidence cutoff
#'
#' Build a line plot of precision, recall, and F1 across confidence cutoffs
#' from a PR-curve data frame (for example, output from
#' \code{compute_variant_pr_curve()}, \code{compute_tree_pr_curve()}, or
#' \code{compute_pr_from_jaccard_mat()}).
#'
#' @param pr_df Data frame containing at least columns
#'   \code{conf_cutoff}, \code{precision}, \code{recall}, and \code{f1}. If
#'   present, \code{fraction_well_matched} is used for optional annotation.
#' @param sample_name Optional title string shown above the panel.
#' @param cutoff Optional numeric cutoff to highlight with a vertical dashed
#'   line and label near the top of the panel.
#' @param x_axis_scale X-axis scale: \code{"log10"} (default) or
#'   \code{"linear"}.
#' @param log10_min Lower bound used for log10 plotting so zero cutoffs can be
#'   displayed safely.
#' @param j_thres Optional Jaccard threshold used only for display annotation
#'   with \code{fraction_well_matched}.
#' @param legend Logical; if \code{TRUE}, show metric legend on the right.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
plot_prec_recall_vs_conf <- function(
    pr_df,
    sample_name = NULL,
    cutoff = NULL,
    x_axis_scale = c("log10", "linear"),
    log10_min = 1e-6,
    j_thres = NULL,
    legend = TRUE
) {
    x_axis_scale <- match.arg(x_axis_scale)
    if (!all(c("conf_cutoff", "precision", "recall", "f1") %in% colnames(pr_df))) {
        stop("pr_df must contain columns: conf_cutoff, precision, recall, f1")
    }

    fraction <- pr_df$fraction_well_matched
    fraction <- if (length(fraction) == 0L) NA_real_ else fraction[1]
    
    if (!is.null(j_thres)) {
        ann <- if (is.na(fraction)) NA_character_ else glue("{signif(100 * fraction, 2)}% vars J >= {j_thres}")
    } else {
        ann = ''
    }

    df_long <- rbind(
        data.frame(conf_cutoff = pr_df$conf_cutoff, metric = "Precision", value = pr_df$precision),
        data.frame(conf_cutoff = pr_df$conf_cutoff, metric = "Recall", value = pr_df$recall),
        data.frame(conf_cutoff = pr_df$conf_cutoff, metric = "F1", value = pr_df$f1)
    )

    x_lab <- "Confidence cutoff"
    vline_x <- cutoff
    x_var <- "conf_cutoff"
    if (x_axis_scale == "log10") {
        df_long$conf_cutoff_plot <- pmax(df_long$conf_cutoff, log10_min)
        x_var <- "conf_cutoff_plot"
        x_lab <- "Confidence cutoff (log10)"
        if (!is.null(vline_x)) vline_x <- max(vline_x, log10_min)
    }

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[x_var]], y = value, color = metric)) +
        ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
        ggplot2::scale_color_manual(
            values = c("Precision" = "#1F77B4", "Recall" = "#D62728", "F1" = "#2CA02C"),
            name = 'Metric'
        ) +
        ggplot2::labs(
            title = sample_name,
            x = x_lab,
            y = NULL
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 10, hjust = 0, margin = margin(b = 10)),
            axis.title = ggplot2::element_text(size = 8),
            axis.text = ggplot2::element_text(size = 7),
            legend.position = if (legend) "right" else "none",
            panel.grid.minor = ggplot2::element_blank()
        )

    if (!is.null(vline_x)) {
        p <- p + ggplot2::geom_vline(
            xintercept = vline_x,
            color = "grey40",
            linetype = "11",
            linewidth = 0.6
        )

        cutoff_label_value <- if (!is.null(cutoff)) cutoff else vline_x
        cutoff_label <- paste0(signif(cutoff_label_value, 3))
        p <- p + ggplot2::annotate(
            "text",
            x = vline_x,
            y = Inf,
            label = cutoff_label,
            vjust = -0.5,
            hjust = 0.5,
            size = 2.8,
            color = "grey25"
        )
    }

    if (!is.na(ann)) {
        x_anno <- min(df_long[[x_var]], na.rm = TRUE)
        y_anno <- min(df_long$value, na.rm = TRUE)
        p <- p + ggplot2::annotate(
            "text",
            x = x_anno,
            y = y_anno,
            label = ann,
            hjust = 0,
            vjust = 0,
            size = 2.6
        )
    }

    if (x_axis_scale == "log10") {
        p <- p + ggplot2::scale_x_continuous(trans = "log10")
    }

    p <- p +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5))

    p
}
