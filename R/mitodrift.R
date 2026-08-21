#' @import dplyr
#' @import tidygraph
#' @import stringr
#' @import ggplot2
#' @import ggtree
#' @importFrom glue glue
#' @importFrom igraph vcount ecount E V V<- E<- 
#' @importFrom phangorn upgma 
#' @importFrom ape root drop.tip nj
#' @importFrom parallelDist parDist
#' @importFrom stats na.omit reorder setNames
#' @useDynLib mitodrift
NULL

#' Optimize tree topology using C++ NNI moves
#'
#' Performs nearest-neighbor interchange (NNI) hill-climbing to find the
#' tree topology that maximizes the belief-propagation score. Uses compiled
#' C++ routines for speed.
#'
#' @param tree_init A rooted `phylo` object used as the starting tree.
#'   Ignored when resuming from an existing trace.
#' @param logP A list of log-probability vectors (one per locus), as returned
#'   by [convert_logliks_to_logP_list()] or [convert_logliks_to_logP_list_colmajor()].
#' @param logA A numeric vector (or list of vectors) of log transition
#'   probabilities, flattened column-major from the transition matrix.
#' @param max_iter Integer; maximum number of NNI iterations.
#' @param outfile Optional file path for saving the tree trace (qs2 format).
#' @param resume Logical; if `TRUE` and `outfile` exists, resume from the
#'   last saved tree instead of starting fresh.
#' @param ncores Integer; number of threads for parallel NNI scoring.
#' @param trace_interval Integer; save the trace to `outfile` every this many
#'   iterations.
#' @return A `multiPhylo` list of trees visited during optimization, each
#'   carrying a `logZ` element with the log-partition-function score.
#' @keywords internal
optimize_tree_cpp = function(
    tree_init = NULL, logP, logA, max_iter = 100,
    outfile = NULL, resume = FALSE, ncores = 1, trace_interval = 5
) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    if (is.list(logA)) {
        score_tree_func = score_tree_bp_wrapper_multi
        nni_func = nni_cpp_parallel_multi
    } else {
        score_tree_func = score_tree_bp_wrapper2
        nni_func = nni_cpp_parallel_cached
    }

    if (resume & !is.null(outfile) & file.exists(outfile)) {
        message('Resuming from existing tree list')
        tree_list = safe_read_chain(outfile)
        tree_current = tree_list %>% .[[length(.)]]
        max_current = tree_current$logZ
        start_iter = length(tree_list)
    } else {
        if (is.null(tree_init) || is.null(logP) || is.null(logA)) {
            stop("tree_init, logP, and logA must be provided when resume is FALSE")
        }
        start_iter = 1
        tree_init = reorder_phylo(tree_init)
        tree_init$edge.length = NULL
        tree_current = tree_init
        max_current = sum(score_tree_func(tree_current$edge, logP, logA))
        tree_current$logZ = max_current
        tree_list = phytools::as.multiPhylo(tree_current)
    }

    runtime = c(0,0,0)

    if (start_iter > max_iter) {
        message(paste0("Already completed ", start_iter - 1, " iterations. No new iterations to run for max_iter = ", max_iter, "."))
        if (!is.null(outfile)) {
            qs2::qd_save(tree_list, outfile)
        }
        return(tree_list)
    }

    for (i in start_iter:max_iter) {

        ptm = proc.time()

        message(paste(i, round(max_current, 4), paste0('(', signif(unname(runtime[3]),2), 's', ')')))

        scores = nni_func(tree_current$edge, logP, logA)
        
        if (max(scores) > max_current) {
            max_id = which.max(scores)
            if (max_id %% 2 == 0) {pair_id = 2} else {pair_id = 1}
            tree_current$edge = matrix(nnin_cpp(tree_current$edge, ceiling(max_id/2))[[pair_id]], ncol = 2)
            tree_current$logZ = max_current = max(scores)
        } else {
            break()
        }

        tree_list = tree_list %>% c(phytools::as.multiPhylo(tree_current))

        if (!is.null(outfile)) {
            if (i == 1 | i %% trace_interval == 0) {
                qs2::qd_save(tree_list, outfile)
            }
        }

        runtime = proc.time() - ptm
        
    }
    
    if (!is.null(outfile)) {
        qs2::qd_save(tree_list, outfile)
    }

    return(tree_list)
}

#' Reorder a phylo object to postorder
#'
#' Creates a deep copy of the phylogeny and reorders its edge matrix to
#' postorder using the compiled C++ helper `reorderRcpp`.
#'
#' @param phy A `phylo` object.
#' @return A new `phylo` object with edges in postorder.
#' @keywords internal
reorder_phylo = function(phy) {
    phy_new = rlang::duplicate(phy, shallow = FALSE)
    phy_new$edge = reorderRcpp(phy$edge) %>% matrix(ncol = 2)
    return(phy_new)
}

#' Convert log-likelihood matrices to a log-probability list (row-major)
#' @keywords internal
#'
#' Like [convert_liks_to_logP_list()] but expects inputs already on the log
#' scale. Produces flat log-probability vectors in row-major layout.
#'
#' @param logliks Named list of log-likelihood matrices (one per variant),
#'   each of dimension `k x n_cells`.
#' @param phy A `phylo` object whose tip labels determine column ordering.
#' @return A named list of numeric vectors, each of length `k * n_nodes`.
convert_logliks_to_logP_list <- function(logliks, phy) {
    
    E <- reorder_phylo(phy)$edge
    phy$node.label <- NULL
    
    P_all <- lapply(logliks, function(liks_mut) {
        n_tips <- length(phy$tip.label)
        n_nodes <- phy$Nnode
        root_node <- E[nrow(E), 1]
        k <- nrow(liks_mut)
        
        P <- matrix(nrow = k, ncol = n_tips + n_nodes)
        rownames(P) <- rownames(liks_mut)

        # Tip likelihoods, internal node likelihoods, root node set up
        P[, 1:n_tips] <- liks_mut[, phy$tip.label]
        P[, (n_tips + 1):(n_tips + n_nodes)] <- log(1/k)
        P[, root_node] <- log(c(1, rep(0, k - 1)))
        
        return(P)
    })

    logP_list <- lapply(P_all, function(P) {
        as.vector(t(P))
    })
    
    return(logP_list)
}

#' Generate VAF bin midpoints
#'
#' Creates `k + 2` evenly spaced VAF bins spanning \[0, 1\] (including
#' boundary bins at 0 and 1) and returns their midpoints.
#'
#' @param k Integer; number of interior VAF bins. The total number of bins
#'   is `k + 2`.
#' @return Numeric vector of bin midpoints of length `k + 2`.
#' @keywords internal
get_vaf_bins = function(k) {
    bins = seq(0, 1, 1/k)
    bins = c(0,bins,1)
    vafs = sapply(1:(length(bins)-1), function(i){(bins[i] + bins[i+1])/2})
    return(vafs)
}

# Caching environments for transition matrices
.mitodrift_T_cache <- new.env(parent = emptyenv())
.mitodrift_Tmat_cache <- new.env(parent = emptyenv())

#' Get the transition matrix for WF model with HMM (with caching)
#' TODO: add log option for small probabilities
#' @param k number of VAF bins
#' @param eps Variant detection error rate
#' @param N population size
#' @param ngen number of generations
#' @param safe whether to add small probability to avoid 0s
#' @return transition matrix
#' @keywords internal
get_transition_mat_wf_hmm <- function(k, eps, N, ngen, safe = FALSE) {
	# Precompute bin boundaries and VAF bin midpoints
	bin_boundaries <- c(0, seq(0, 1, 1 / k), 1)
	vaf_bins <- get_vaf_bins(k)

	# Zero generations: identity (then boundary rows adjusted below)
	if (ngen == 0) {
		A <- diag(k + 2)
	} else {
		# ---- T_ngen cache key depends only on N (matrix size) and ngen (power) ----
		key_T <- paste0(N, "|", as.integer(ngen))

		if (exists(key_T, envir = .mitodrift_T_cache, inherits = FALSE)) {
			T_ngen <- get(key_T, envir = .mitodrift_T_cache, inherits = FALSE)
		} else {
			# Build or reuse the single-generation transition matrix T_mat for this N
			key_Tmat <- as.character(N)
			if (exists(key_Tmat, envir = .mitodrift_Tmat_cache, inherits = FALSE)) {
				T_mat <- get(key_Tmat, envir = .mitodrift_Tmat_cache, inherits = FALSE)
			} else {
				p_vec <- (0:N) / N
				T_mat <- outer(p_vec, 0:N, function(p, k) dbinom(k, size = N, prob = p))
				assign(key_Tmat, T_mat, envir = .mitodrift_Tmat_cache)
			}

			# Multi-generation transition via integer matrix power
			T_ngen <- expm::`%^%`(T_mat, ngen)
			assign(key_T, T_ngen, envir = .mitodrift_T_cache)
		}

		# ---- Aggregate allele-count probabilities into VAF bins ----
		A <- matrix(0, nrow = k + 2, ncol = k + 2)
		for (i in 1:(k + 2)) {
			# Representative starting allele count for bin i (midpoint mapping)
			start_vaf <- vaf_bins[i]
			start_allele_count <- round(start_vaf * N)
			start_allele_count <- max(0, min(N, start_allele_count))

			# Probability distribution after ngen generations
			prob_dist_after_ngen <- T_ngen[start_allele_count + 1, ]

			for (j in 1:(k + 2)) {
				xstart <- floor(bin_boundaries[j] * N) + 1
				xend <- floor(bin_boundaries[j + 1] * N)
				xstart <- min(xstart, xend)
				if (j == k + 1) xend <- xend - 1
				if (xstart <= xend) {
					indices <- xstart:xend
					A[i, j] <- sum(prob_dist_after_ngen[indices + 1])
				} else {
					A[i, j] <- 0
				}
			}
		}
	}

	# ---- Apply mutation rate to boundary rows ----
	A[1, ] <- c(1 - eps, rep(eps / (k + 1), k + 1))
	A[nrow(A), ] <- rev(A[1, ])

	# Names and safety floor
	colnames(A) <- vaf_bins
	rownames(A) <- vaf_bins
	if (safe) {
		A[A == 0] <- 1e-16
	}
	return(A)
}

#' Wrapper function to interpolate non-integer generations
#' @param k number of VAF bins
#' @param eps Variant detection error rate
#' @param N population size
#' @param ngen number of generations (may be non-integer)
#' @param safe Logical; if `TRUE`, replace zero entries with a small floor
#'   value to avoid numerical issues.
#' @return transition matrix
#' @keywords internal
#' @noRd
get_transition_mat_wf_hmm_wrapper = function(k, eps, N, ngen, safe = FALSE) {
    
    # Check if ngen is an integer
    if (ngen == round(ngen)) {
        # If integer, use the original function directly
        return(get_transition_mat_wf_hmm(k = k, eps = eps, N = N, ngen = ngen, safe = safe))
    }
    
    # If non-integer, interpolate between two nearest integer generations
    ngen_lower = floor(ngen)
    ngen_upper = ceiling(ngen)
    
    # Get transition matrices for the two nearest integer generations
    A_lower = get_transition_mat_wf_hmm(k = k, eps = eps, N = N, ngen = ngen_lower, safe = safe)
    A_upper = get_transition_mat_wf_hmm(k = k, eps = eps, N = N, ngen = ngen_upper, safe = safe)
    
    # Calculate interpolation weight
    weight = ngen - ngen_lower
    
    # Linear interpolation between the two matrices
    A_interpolated = (1 - weight) * A_lower + weight * A_upper
    
    return(A_interpolated)
}

#' Make a rooted NJ tree
#' @param vmat A matrix of cell-by-variable values
#' @param dist_method The distance method to use
#' @param ncores Number of threads for `parallelDist::parDist` (default: 1)
#' @return A phylo object
#' @keywords internal
make_rooted_nj = function(vmat, dist_method = 'manhattan', ncores = 1) {
    vmat[is.na(vmat)] = 0
    vmat = cbind(vmat, outgroup = 0) %>% as.matrix %>% t
    if (ncores > 1) {
        dist_mat = parallelDist::parDist(vmat, method = dist_method, threads = ncores)
    } else {
        dist_mat = dist(vmat, method = dist_method)
    }
    nj_tree = ape::nj(dist_mat) %>% 
        ape::root(outgroup = 'outgroup') %>% drop.tip('outgroup') 
    return(nj_tree)
}

#' Decode a tree using CRF belief propagation (R version)
#'
#' Constructs a conditional random field (CRF) on the tree with a single
#' shared transition matrix and computes per-variant marginal beliefs via
#' tree belief propagation. Optionally returns posterior means, MAP
#' assignments, or the full CRF objects.
#'
#' @param tn A `phylo` object representing the tree topology.
#' @param A Transition matrix (square, `k x k`) with VAF bin midpoints as
#'   row/column names.
#' @param liks Named list of likelihood matrices (one per variant), each
#'   `k x n_cells`.
#' @param post_max Logical; if `TRUE`, also compute MAP (Viterbi) decoding.
#' @param store_bels Logical; if `TRUE`, store per-variant node and edge
#'   beliefs in the output.
#' @param store_crfs Logical; if `TRUE`, store a copy of the CRF object for
#'   each variant.
#' @param debug Logical; if `TRUE`, return a detailed list instead of just
#'   the `tbl_graph`.
#' @param score_only Logical; if `TRUE`, skip posterior computation and only
#'   attach log-partition scores.
#' @return A `tbl_graph` tree with per-variant posterior means (columns
#'   `p_<variant>`) and a `logZ` vector of log-partition-function values.
#'   When `debug = TRUE`, a list with additional diagnostic components.
#' @keywords internal
decode_tree = function(
    tn, A, liks, post_max = FALSE, store_bels = FALSE, store_crfs = FALSE, debug = FALSE,
    score_only = FALSE
) {
    
    k = ncol(A)
    vafs = as.numeric(colnames(A))
    # convert tree to CRF
    tn$node.label = NULL
    Gn = as.igraph(tn)
    
    gtree = as_tbl_graph(Gn)
    root_node = gtree %>% filter(node_is_root()) %>% pull(name) %>% as.character

    adj_n = as_adjacency_matrix(Gn)
    crf = make.crf(adj_n, k)

    # add edge potentials
    flip_dict = Gn %>% 
        as_edgelist(names = F) %>% 
        apply(1, function(x){
            setNames(is.unsorted(x), paste0(sort(x), collapse = ','))
        }, simplify = F) %>%
        unlist

    for (i in 1:crf$n.edges) {
        
        epair = paste0(sort(crf$edges[i,]), collapse = ',')
        flip = flip_dict[[epair]]
        
        if (flip) {
            crf$edge.pot[[i]] = t(A)
        } else {
            crf$edge.pot[[i]] = A
        }
        
    }

    # add note potentials
    vnames = names(V(Gn))
    crf$node.labels = vnames
    rownames(crf$node.pot) = vnames

    logZ = c()
    ebels = list()
    nbels = list()
    crfs = list()
    res_max = NULL

    for (mut in names(liks)) {

        crf$node.pot[colnames(liks[[mut]]),] = t(liks[[mut]])
        crf$node.pot[!vnames %in% colnames(liks[[mut]]),] = 1/k
        crf$node.pot[root_node,] = c(1, rep(0, k-1))
        
        # decoding
        res_mar = infer.tree(crf)
        
        if (!score_only) {
            # append posterior mean to graph tree
            p_dat = res_mar$node.bel %*% diag(vafs) %>% rowSums
            p_dat = p_dat %>% data.frame(vnames, .) %>%
                setNames(c('name', paste0('p_', mut)))
            
            gtree = gtree %>% activate(nodes) %>%
                select(-any_of(c(paste0('p_', mut)))) %>% 
                left_join(p_dat, by = join_by(name))

            if (post_max) {
                res_max = decode.tree(crf)
                z_dat = data.frame(vnames, vafs[res_max]) %>% setNames(c('name', paste0('z_', mut))) 

                gtree = gtree %>% activate(nodes) %>%
                    select(-any_of(c(paste0('z_', mut)))) %>% 
                    left_join(z_dat, by = join_by(name))
            }
        }

        logZ = c(logZ, res_mar$logZ)

        if (store_bels) {
            ebels[[mut]] = res_mar$edge.bel
            nbels[[mut]] = res_mar$node.bel
            rownames(nbels[[mut]]) = vnames
        }

        if (store_crfs) {
            crf_copy <- rlang::env_clone(crf)
            attributes(crf_copy) <- attributes(crf)
            crfs[[mut]] = crf_copy
        }
    }

    logZ = setNames(logZ, names(liks))
    gtree$logZ = logZ

    if (debug) {
        return(list('gtree' = gtree, 'crfs' = crfs, 'Gn' = Gn, 
        'res_mar' = res_mar, 'res_max' = res_max, 'ebels' = ebels, 'nbels' = nbels))
    }

    return(gtree)
}

################################### MCMC ######################################

#' Attach a new edge matrix to a phylo object
#'
#' Replaces the edge matrix of a `phylo` object with a new one (e.g. from
#' an MCMC sample).
#'
#' @param phy A `phylo` object serving as the template.
#' @param edges Integer vector to be reshaped into a 2-column edge matrix.
#' @return A `phylo` object with the updated edge matrix.
#' @keywords internal
#' @noRd
attach_edges = function(phy, edges) {

    phy_new = phy
    E_new = matrix(edges, ncol = 2)
    phy_new$edge = E_new

    return(phy_new)
}

#' Safely read a qs2 chain file
#'
#' Reads a chain trace with error handling for truncated or missing files.
#' Understands two on-disk layouts and combines them transparently so every
#' caller gets the same list-of-chains-of-trees object regardless of which
#' one produced it:
#'   - legacy: a single qs2 file at `path` holding the complete trace
#'     (including the seed tree as each chain's first element). This is the
#'     only format archived regression keys (`test_dat/keys/`) ever use, and
#'     it is never rewritten or migrated -- only read.
#'   - block-wise: a `<path>.blocks/` directory of small, numbered qs2 files,
#'     each holding only the iterations added at that checkpoint. Written by
#'     `run_tree_mcmc_batch()` going forward so each checkpoint is an append
#'     (cost proportional to that checkpoint's batch, not the whole trace).
#' If both exist (e.g. resuming a run that was started before block-wise
#' writing existed), the legacy file is treated as an implicit first block
#' and the numbered blocks are appended after it in order.
#'
#' @param path Character file path (same `outfile`/`trace_file` path used
#'   everywhere else; the block directory, if any, is derived from it).
#' @param ncores Integer; number of threads for `qs2::qd_read`.
#' @return A list of per-chain edge-list traces, or `NULL` if nothing is
#'   found at `path` in either layout.
#' @keywords internal
#' @noRd
trace_blocks_dir = function(path) paste0(path, '.blocks')

#' @keywords internal
#' @noRd
list_trace_blocks = function(path) {
    blocks_dir = trace_blocks_dir(path)
    if (!dir.exists(blocks_dir)) return(character(0))
    files = list.files(blocks_dir, pattern = '^[0-9]{8}\\.qs2$', full.names = TRUE)
    sort(files)
}

#' @keywords internal
#' @noRd
read_qs2_file = function(path, ncores = 1) {
    if (!file.exists(path) || dir.exists(path)) return(NULL)
    fi = file.info(path)
    if (is.na(fi$size) || fi$size <= 0) return(NULL)
    tryCatch(qs2::qd_read(path, nthreads = ncores), error = function(e) NULL)
}

#' @keywords internal
#' @noRd
safe_read_chain = function(path, ncores = 1) {
    legacy = read_qs2_file(path, ncores = ncores)
    block_files = list_trace_blocks(path)
    if (length(block_files) == 0) return(legacy)

    blocks = lapply(block_files, read_qs2_file, ncores = ncores)
    blocks = Filter(Negate(is.null), blocks)
    parts = if (is.null(legacy)) blocks else c(list(legacy), blocks)
    if (length(parts) == 0) return(NULL)

    nchains_here = max(vapply(parts, length, integer(1)))
    combined = vector('list', nchains_here)
    for (chain_id in seq_len(nchains_here)) {
        combined[[chain_id]] = do.call(c, lapply(parts, function(part) {
            if (chain_id <= length(part) && !is.null(part[[chain_id]])) part[[chain_id]] else list()
        }))
    }
    names(combined) = as.character(seq_len(nchains_here))
    combined
}

#' Run tree-topology MCMC in batches with convergence monitoring
#'
#' Runs MCMC sampling in fixed-size batches, computing ASDSF convergence
#' diagnostics between batches. Supports automatic stopping when ASDSF
#' drops below a threshold and resume from a previous run.
#'
#' @param phy_init A rooted `phylo` object used as the starting tree.
#' @param logP_list List of log-probability vectors (one per locus).
#' @param logA_vec Numeric vector of log transition probabilities.
#' @param outfile File path identifying the trace (qs2 format). Checkpoints
#'   are written as small per-checkpoint block files under `<outfile>.blocks/`
#'   rather than rewriting one cumulative file; a legacy single file directly
#'   at `outfile` (as archived regression keys use) is also still read
#'   transparently. Use `safe_read_chain()` rather than `qs2::qd_read()`
#'   directly to read a trace written by this function.
#' @param diagfile Optional file path for saving convergence diagnostics
#'   (RDS format).
#' @param diag Logical; whether to compute diagnostics (currently unused,
#'   diagnostics are always computed).
#' @param max_iter Integer; total number of MCMC iterations per chain
#'   (ignored when `conv_thres` is set).
#' @param nchains Integer; number of independent chains.
#' @param ncores Integer; number of threads for C++ MCMC sampling.
#' @param ncores_qs Integer; number of threads for qs2 serialization.
#' @param batch_size Integer; number of iterations per batch.
#' @param conv_thres Numeric or `NULL`; if set, run until the ASDSF drops
#'   below this threshold instead of using `max_iter`.
#' @param resume Logical; if `TRUE`, resume from existing `outfile`.
#' @param mc3_temperatures Optional numeric temperature ladder for
#'   Metropolis-coupled MCMC. It must start at 1 and strictly increase. `NULL`
#'   or the single value `1` runs the original independent-chain sampler.
#' @param mc3_swap_interval Positive integer; iterations between proposed
#'   swaps of adjacent temperatures.
#' @param mc3_swap_scheme MC3 adjacent-swap schedule. `"rnn"` preserves the
#'   original randomly selected adjacent-pair proposal at each swap barrier;
#'   `"deo"` alternates all disjoint even and odd adjacent pairs.
#' @param mc3_ncores Optional number of sampling threads for MC3. Defaults to
#'   `ncores`; set it up to `nchains * length(mc3_temperatures)` to update all
#'   ensemble-temperature pairs concurrently when those CPUs are allocated.
#' @param mc3_statefile Optional RDS checkpoint for heated-chain states and
#'   swap statistics. Defaults to `paste0(outfile, ".mc3_state.rds")`. Also
#'   carries an `asdsf_state` field (the incremental target-clade ASDSF
#'   accumulator's tiny `counts`/`totals` state, not the trace) once a run has
#'   completed at least one batch -- when present on an MC3 resume, this lets
#'   the run skip reading the full historical trace back into memory
#'   entirely, since MC3 continuation state lives here, not in the cold-chain
#'   trace, so nothing from the trace itself is otherwise needed at resume
#'   time. Checkpoints from before this field existed simply lack it, so
#'   resume falls back to reconstructing ASDSF state from the full trace, as
#'   before -- fully backward compatible.
#' @param checkpoint_every Positive integer; persist a new checkpoint and the
#'   sampler state every this many diagnostic batches. The final batch is always
#'   persisted. Values greater than one reduce I/O at the cost of re-running up
#'   to `checkpoint_every - 1` batches after interruption. Each checkpoint is
#'   written as its own small block file (see `outfile`), so this cost no
#'   longer grows with total run length.
#' @param return_trace Logical; if `FALSE`, skip reconstructing the full
#'   trace for the return value and return `NULL` instead. The checkpointed
#'   file(s) on disk are unaffected either way. Set `FALSE` when the caller
#'   only needs the on-disk trace (e.g. to avoid paying for a read-back that
#'   would just be discarded).
#' @return A list of edge-list chains (one list of edge matrices per chain),
#'   or `NULL` when `return_trace = FALSE`.
#' @keywords internal
run_tree_mcmc_batch = function(
    phy_init, logP_list, logA_vec, outfile, diagfile = NULL, diag = TRUE, max_iter = 10000, nchains = 1, ncores = 1, ncores_qs = 1,
    batch_size = 1000, conv_thres = NULL, resume = FALSE,
    mc3_temperatures = NULL, mc3_swap_interval = 10L,
    mc3_swap_scheme = c('rnn', 'deo'), mc3_ncores = NULL,
    mc3_statefile = NULL,
    checkpoint_every = 1L,
    return_trace = TRUE
) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)

    ncores_qs <- if (isTRUE(qs2:::check_TBB())) ncores_qs else 1L
    message('Using ', ncores_qs, ' cores for saving/writing MCMC trace')
    checkpoint_every <- as.integer(checkpoint_every)
    if (length(checkpoint_every) != 1L || is.na(checkpoint_every) || checkpoint_every < 1L) {
        stop('checkpoint_every must be a positive integer')
    }

    mc3_swap_scheme <- match.arg(mc3_swap_scheme)
    mc3_requested <- !is.null(mc3_temperatures)
    use_mc3 <- FALSE
    if (mc3_requested) {
        mc3_temperatures <- as.numeric(mc3_temperatures)
        if (length(mc3_temperatures) < 1L || mc3_temperatures[[1]] != 1 ||
            any(!is.finite(mc3_temperatures)) || any(diff(mc3_temperatures) <= 0)) {
            stop('mc3_temperatures must contain finite, strictly increasing values starting at 1')
        }
        use_mc3 <- length(mc3_temperatures) > 1L
        mc3_swap_interval <- as.integer(mc3_swap_interval)
        if (length(mc3_swap_interval) != 1L || is.na(mc3_swap_interval) || mc3_swap_interval < 1L) {
            stop('mc3_swap_interval must be a positive integer')
        }
        if (use_mc3) {
            if (is.null(mc3_statefile)) mc3_statefile <- paste0(outfile, '.mc3_state.rds')
        } else {
            message('A single MC3 temperature (1) uses the ordinary MCMC sampler')
        }
    }
    sampling_ncores <- if (use_mc3 && !is.null(mc3_ncores)) mc3_ncores else ncores
    sampling_ncores <- as.integer(sampling_ncores)
    if (length(sampling_ncores) != 1L || is.na(sampling_ncores) || sampling_ncores < 1L) {
        stop('mc3_ncores must be a positive integer when supplied')
    }
    RcppParallel::setThreadOptions(numThreads = sampling_ncores)
    if (use_mc3) {
        message('MC3 enabled with temperatures: ', paste(mc3_temperatures, collapse = ', '),
                '; adjacent swap interval: ', mc3_swap_interval,
                '; swap scheme: ', mc3_swap_scheme,
                '; sampling threads: ', sampling_ncores)
    }

    if (resume && !is.null(diagfile) && file.exists(diagfile)) {
        diag_history = readRDS(diagfile)
    } else {
        diag_history = data.frame()
    }

    chains = 1:nchains

    outdir = dirname(outfile)

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    blocks_dir = trace_blocks_dir(outfile)
    if (!resume) {
        # A fresh run must never build on whatever happens to already be at
        # `outfile` -- e.g. a leftover trace from an earlier, unrelated run
        # at the same path. Clear both possible layouts so this run starts
        # from a genuinely empty slate rather than silently combining with
        # stale data (the block-wise write path below never reads existing
        # content to decide what to write, so this check is the only thing
        # standing between a fresh run and exactly that contamination).
        if (file.exists(outfile) && !dir.exists(outfile)) file.remove(outfile)
        if (dir.exists(blocks_dir)) unlink(blocks_dir, recursive = TRUE)
    }

    # MC3 fast-resume: if the heated-state checkpoint already carries an asdsf_state
    # field (written below on any run that completed at least one batch under this
    # code), skip the full-trace read entirely. MC3 continuation state lives entirely
    # in mc3_checkpoint$states (used directly by tree_mc3_parallel_seeded() below); the
    # cold-chain trace itself is not consulted by the MC3 sampling path at all (unlike
    # the ordinary-MCMC path, which does need the trace's last tree per chain as its
    # starting point). Checkpoints from before this field existed simply lack it, so
    # this falls back to the slow path below with no special-case handling needed.
    mc3_checkpoint <- NULL
    asdsf_fastpath <- FALSE
    if (use_mc3 && resume && file.exists(mc3_statefile)) {
        mc3_checkpoint <- readRDS(mc3_statefile)
        checkpoint_swap_scheme <- if (is.null(mc3_checkpoint$swap_scheme)) {
            'rnn'
        } else {
            as.character(mc3_checkpoint$swap_scheme)
        }
        if (!identical(as.numeric(mc3_checkpoint$temperatures), mc3_temperatures) ||
            !identical(as.integer(mc3_checkpoint$swap_interval), mc3_swap_interval) ||
            !identical(checkpoint_swap_scheme, mc3_swap_scheme)) {
            stop('MC3 checkpoint temperature ladder, swap interval, or swap scheme does not match the requested run')
        }
        # Checkpoints written before swap schemes were configurable used RNN.
        # Persist that explicit interpretation at the next checkpoint write.
        mc3_checkpoint$swap_scheme <- checkpoint_swap_scheme
        asdsf_fastpath <- !is.null(mc3_checkpoint$asdsf_state)
    }

    if (asdsf_fastpath) {
        message('MC3 fast-resume: asdsf state found in checkpoint, skipping full-trace read ',
                '(peak memory no longer grows with total accumulated iterations on resume)')
        edge_list_all <- NULL  # deliberately not loaded -- last_tree (derived from it) is
                                # unused on the MC3 sampling path, see below
        completed_iters <- as.integer(mc3_checkpoint$completed_iters)
        asdsf_state <- mc3_checkpoint$asdsf_state
    } else {
        edge_list_all = if (resume) safe_read_chain(outfile, ncores = ncores_qs) else NULL
        if (is.null(edge_list_all)) {
            edge_list_all = vector('list', nchains)
        } else {
            length(edge_list_all) = nchains
        }

        max_len = if (is.null(conv_thres)) max_iter + 1L else NULL
        for (i in seq_along(edge_list_all)) {
            chain_list = edge_list_all[[i]]
            if (is.null(chain_list) || length(chain_list) == 0) {
                edge_list_all[[i]] = list(phy_init$edge)
                next
            }
            if (!is.null(max_len) && length(chain_list) > max_len) {
                edge_list_all[[i]] = chain_list[seq_len(max_len)]
            }
        }
        names(edge_list_all) = as.character(chains)

        if (!is.null(mc3_checkpoint)) {
            checkpoint_iter <- if (is.null(mc3_checkpoint$completed_iters)) NA_integer_ else as.integer(mc3_checkpoint$completed_iters)
            trace_iter <- length(edge_list_all[[1]]) - 1L
            if (is.na(checkpoint_iter) || checkpoint_iter != trace_iter) {
                stop('MC3 heated-state checkpoint is not synchronized with the cold-chain trace')
            }
        }

        completed_iters = length(edge_list_all[[1]]) - 1L
        asdsf_state <- initialize_target_tree_asdsf_state(
            phy_init, edge_list_all, rooted = TRUE)
    }
    if (use_mc3 && is.null(mc3_checkpoint)) {
        if (resume && any(lengths(edge_list_all) > 1L)) {
            stop('Cannot resume MC3: heated-chain checkpoint is missing: ', mc3_statefile)
        }
        mc3_checkpoint <- list(
            states = lapply(chains, function(i) rep(list(phy_init$edge), length(mc3_temperatures))),
            temperatures = mc3_temperatures,
            swap_interval = mc3_swap_interval,
            swap_scheme = mc3_swap_scheme,
            deo_phase = 0L,
            swap_attempts = matrix(0L, nrow = nchains, ncol = length(mc3_temperatures) - 1L),
            swap_accepts = matrix(0L, nrow = nchains, ncol = length(mc3_temperatures) - 1L),
            completed_iters = 0L
        )
    }

    if (nrow(diag_history) > 0L && 'completed_iters' %in% names(diag_history)) {
        diag_history <- diag_history[diag_history$completed_iters <= completed_iters, , drop = FALSE]
    }

    # Early-return trace value: reuse the in-memory edge_list_all when the slow path
    # already loaded it (byte-identical to what it always returned); the fast path
    # doesn't have it in memory, so read it fresh here rather than returning NULL --
    # this is the rare instant-already-done case, not the steady-state resume this
    # optimization targets, so paying for one read here doesn't defeat the point.
    get_return_trace <- function() {
        if (!is.null(edge_list_all)) edge_list_all else safe_read_chain(outfile, ncores = ncores_qs)
    }

    if (is.null(conv_thres)) {
        remaining = max_iter - completed_iters
        message('Remaining iterations per chain: ', remaining)
        if (remaining <= 0L) {
            message('All chains have completed the requested iterations.')
            return(if (return_trace) get_return_trace() else invisible(NULL))
        }
        total_batches = ceiling(remaining / batch_size)
        message('Running MCMC with ', length(chains), ' chains in up to ', total_batches, ' batches of ', batch_size)
    } else {
        if (nrow(diag_history) > 0) {
            last_asdsf <- utils::tail(diag_history$asdsf, 1)
            if (last_asdsf <= conv_thres) {
                message('Convergence threshold (ASDSF) already reached (', signif(last_asdsf, 4), '). Nothing to run.')
                return(if (return_trace) get_return_trace() else invisible(NULL))
            } else {
                message('Last recorded ASDSF: ', signif(last_asdsf, 4))
            }
        }
        message('Running MCMC with ', length(chains), ' chains until ASDSF <= ', conv_thres,
                ' (batch size ', batch_size, '; no iteration cap)')
    }

    # From here on the full historical trace is only needed on disk, not in R
    # memory: ASDSF is already updated incrementally per batch (asdsf_state
    # above), the on-disk file is guaranteed up to date whenever the loop
    # exits (persist_batch below always fires on the final batch), and the
    # return value is reconstructed with a single read at the very end.
    # Keeping `edge_list_all` growing for the whole run is unnecessary and,
    # for large trees run to convergence, can end up dwarfing the NNICache's
    # own (fixed-size) memory footprint -- see NNI_CACHE_BIDIRECTIONAL_NOTES.md
    # for the sizing comparison that motivated this change. `last_tree` seeds
    # the next batch; `pending` accumulates only what hasn't been flushed to
    # disk yet, so peak memory is bounded by `checkpoint_every` batches'
    # worth of trace, not by the total iteration count.
    # last_tree seeds tree_mcmc_parallel_seeded()'s starting point on the ordinary-MCMC
    # path; MC3 continuation instead comes entirely from mc3_checkpoint$states, so
    # last_tree is unused whenever use_mc3 is TRUE (in particular, on the fast-resume
    # path above, where edge_list_all is deliberately not loaded and last_tree is simply
    # a placeholder that gets overwritten per-chain in the batch loop regardless).
    last_tree = if (is.null(edge_list_all)) {
        vector('list', nchains)
    } else {
        lapply(edge_list_all, function(chain_list) chain_list[[length(chain_list)]])
    }
    pending = vector('list', nchains)
    for (i in chains) pending[[i]] = list()
    names(pending) = as.character(chains)
    rm(edge_list_all)
    # Block-wise checkpoint state. next_block_idx continues from whatever is
    # already on disk (a real resume) or starts at 1 (the clean slate a
    # fresh run just established above). seed_written tracks whether the
    # seed tree (phy_init$edge) has already been persisted as the first
    # element of some chain's trace -- true on any resume with existing
    # content, false only when starting completely fresh, so the very first
    # block written below includes it and every later block doesn't
    # duplicate it.
    existing_blocks = list_trace_blocks(outfile)
    next_block_idx = length(existing_blocks) + 1L
    seed_written = resume && (file.exists(outfile) || length(existing_blocks) > 0L)

    batch_idx = 0L
    repeat {
        if (is.null(conv_thres)) {
            remaining = max_iter - completed_iters
            if (remaining <= 0L) {
                message('All chains have completed the requested iterations.')
                break
            }
            iter_this_batch = min(batch_size, remaining)
            batch_label = paste('batch', batch_idx + 1L, '; remaining:', ceiling(remaining / batch_size))
        } else {
            iter_this_batch = batch_size
            batch_label = paste('batch', batch_idx + 1L)
        }

        batch_idx = batch_idx + 1L
        message('Running ', batch_label)

        ptm <- proc.time()

        iter_vec = rep(iter_this_batch, length(chains))
        start_edges = last_tree
        seed_index <- if (use_mc3) completed_iters else batch_idx - 1L
        seed_vec = as.integer((1000003 * seed_index + chains) %% .Machine$integer.max)

        if (use_mc3) {
            if (mc3_swap_scheme == 'deo') {
                deo_phase <- if (is.null(mc3_checkpoint$deo_phase)) {
                    0L
                } else {
                    as.integer(mc3_checkpoint$deo_phase)
                }
                mc3_result <- tree_mc3_parallel_seeded_deo(
                    mc3_checkpoint$states,
                    logP_list,
                    logA_vec,
                    iter_vec,
                    seed_vec,
                    mc3_temperatures,
                    mc3_swap_interval,
                    deo_phase
                )
                mc3_checkpoint$deo_phase <- as.integer(
                    (deo_phase + iter_this_batch %/% mc3_swap_interval) %% 2L)
            } else {
                mc3_result <- tree_mc3_parallel_seeded(
                    mc3_checkpoint$states,
                    logP_list,
                    logA_vec,
                    iter_vec,
                    seed_vec,
                    mc3_temperatures,
                    mc3_swap_interval
                )
            }
            elist_active <- mc3_result$traces
            # Keep this batch's raw swap counts separate from the cumulative
            # checkpoint totals.  The latter are required to resume the
            # sampler, while the former make recent/windowed acceptance
            # diagnostics possible without losing the adjacent-pair detail.
            # Aggregate over independent ensembles here: every column is one
            # adjacent temperature pair and the exact integer numerators and
            # denominators remain available in `diag_history` below.
            mc3_batch_swap_attempts <- colSums(mc3_result$swap_attempts)
            mc3_batch_swap_accepts <- colSums(mc3_result$swap_accepts)
            mc3_pair_temperatures <- cbind(
                lower = head(mc3_temperatures, -1L),
                upper = tail(mc3_temperatures, -1L)
            )
            pair_names <- paste0('pair_', seq_along(mc3_batch_swap_attempts))
            names(mc3_batch_swap_attempts) <- pair_names
            names(mc3_batch_swap_accepts) <- pair_names
            rownames(mc3_pair_temperatures) <- pair_names
            mc3_batch_swap_acceptance_by_pair <- ifelse(
                mc3_batch_swap_attempts > 0L,
                mc3_batch_swap_accepts / mc3_batch_swap_attempts,
                NA_real_
            )
            mc3_checkpoint$states <- lapply(mc3_result$final_states, function(ensemble) {
                lapply(ensemble, function(edges) matrix(edges, ncol = 2))
            })
            mc3_checkpoint$swap_attempts <- mc3_checkpoint$swap_attempts + mc3_result$swap_attempts
            mc3_checkpoint$swap_accepts <- mc3_checkpoint$swap_accepts + mc3_result$swap_accepts
        } else {
            elist_active = tree_mcmc_parallel_seeded(
                start_edges,
                logP_list,
                logA_vec,
                iter_vec,
                seed_vec
            )
        }

        new_edge_list_chains <- vector('list', nchains)
        for (chain_id in chains) {
            elist = restore_elist(elist_active[[chain_id]])
            if (length(elist) > 0) {
                elist = elist[-1]
            }
            new_edge_list_chains[[chain_id]] <- elist
            pending[[chain_id]] = c(pending[[chain_id]], elist)
            if (length(elist) > 0) {
                last_tree[[chain_id]] = elist[[length(elist)]]
            }
        }
        completed_iters = completed_iters + iter_this_batch

        asdsf_state <- update_target_tree_asdsf_state(
            asdsf_state, phy_init, new_edge_list_chains, rooted = TRUE)
        asdsf <- target_tree_asdsf_from_state(asdsf_state, min_freq = 0)
        message('ASDSF (target clades) after ', batch_label, ': ', signif(asdsf, 4))
        if (!is.null(diagfile)) {
            diag_entry = data.frame(
                batch = batch_idx,
                completed_iters = completed_iters,
                asdsf = asdsf,
                mc3_swap_acceptance = if (use_mc3) {
                    sum(mc3_checkpoint$swap_accepts) / sum(mc3_checkpoint$swap_attempts)
                } else NA_real_
            )
            if (use_mc3) {
                batch_attempts <- sum(mc3_batch_swap_attempts)
                batch_accepts <- sum(mc3_batch_swap_accepts)
                diag_entry$mc3_swap_attempts_batch <- batch_attempts
                diag_entry$mc3_swap_accepts_batch <- batch_accepts
                diag_entry$mc3_swap_acceptance_batch <- if (batch_attempts > 0L) {
                    batch_accepts / batch_attempts
                } else NA_real_
                diag_entry$mc3_swap_attempts_by_pair_batch <- I(list(mc3_batch_swap_attempts))
                diag_entry$mc3_swap_accepts_by_pair_batch <- I(list(mc3_batch_swap_accepts))
                diag_entry$mc3_swap_acceptance_by_pair_batch <- I(list(mc3_batch_swap_acceptance_by_pair))
                diag_entry$mc3_swap_pair_temperatures <- I(list(mc3_pair_temperatures))
            }
            diag_history = bind_rows(diag_history, diag_entry)
        }
        converged <- !is.null(conv_thres) && !is.na(asdsf) && asdsf <= conv_thres
        fixed_complete <- is.null(conv_thres) && completed_iters >= max_iter
        persist_batch <- batch_idx %% checkpoint_every == 0L || converged || fixed_complete
        if (persist_batch) {
            # Write only what's new since the last checkpoint as its own
            # block -- no read of prior data, so cost is proportional to
            # this checkpoint's batch, not the whole trace so far. The first
            # block of a fresh run also carries the seed tree, matching what
            # a single element of a legacy trace's first position held.
            block_payload = vector('list', nchains)
            for (chain_id in chains) {
                block_payload[[chain_id]] = if (!seed_written) {
                    c(list(phy_init$edge), pending[[chain_id]])
                } else {
                    pending[[chain_id]]
                }
            }
            names(block_payload) = as.character(chains)
            dir.create(blocks_dir, recursive = TRUE, showWarnings = FALSE)
            block_file = file.path(blocks_dir, sprintf('%08d.qs2', next_block_idx))
            tmp_file = paste0(block_file, '.tmp')
            qs2::qd_save(block_payload, tmp_file, nthreads = ncores_qs)
            file.rename(tmp_file, block_file)
            rm(block_payload)
            next_block_idx = next_block_idx + 1L
            seed_written = TRUE
            for (i in chains) pending[[i]] = list()

            if (use_mc3) {
                mc3_checkpoint$completed_iters <- completed_iters
                # asdsf_state is tiny (nchains x n_target_clades counts + a length-nchains
                # totals vector) -- carrying it here is what lets a later resume skip
                # reconstructing it from the full historical trace (see asdsf_fastpath
                # above). No separate file: this is already written/read at exactly the
                # cadence needed, so a second checkpoint file would just be more state to
                # keep in sync for no benefit.
                mc3_checkpoint$asdsf_state <- asdsf_state
                saveRDS(mc3_checkpoint, mc3_statefile)
            }
            if (!is.null(diagfile)) saveRDS(diag_history, diagfile)
        }
        if (converged) {
            message('Convergence threshold (ASDSF) reached. Stopping MCMC.')
            break
        }
        batch_time <- proc.time() - ptm
        message(paste('Completed', batch_label, paste0('(', signif(batch_time[['elapsed']], 2), 's', ')')))
    }

    # persist_batch always fires on the loop's final batch (converged or
    # fixed_complete), so the file is guaranteed to hold the complete,
    # up-to-date trace here -- reconstruct the return value from it rather
    # than from an in-memory accumulator. Skippable via return_trace = FALSE
    # for callers (e.g. Mitodrift$run_mcmc()) that only need the checkpointed
    # file on disk and would otherwise pay for reading the whole trace back
    # just to discard it.
    if (!return_trace) return(invisible(NULL))
    return(safe_read_chain(outfile, ncores = ncores_qs))
}

#' Restore edge-list vectors to 2-column matrices
#'
#' Converts a list of flat integer vectors (from C++ output) back into
#' proper 2-column edge matrices.
#'
#' @param elist List of integer vectors, each representing a flattened
#'   edge matrix.
#' @return List of 2-column integer matrices.
#' @keywords internal
#' @noRd
restore_elist = function(elist) {
    lapply(elist, function(edges){matrix(edges, ncol = 2)})
}

#' Collect MCMC chains into a multiPhylo object
#'
#' Reconstructs `phylo` trees from raw edge-list chains, applies burn-in
#' removal and iteration truncation, then pools all chains into a single
#' `multiPhylo` object.
#'
#' @param edge_list_all List of chains, each a list of 2-column integer
#'   edge matrices.
#' @param phy_init The initial `phylo` object whose tip labels and metadata
#'   are used to reconstruct full `phylo` objects.
#' @param burnin Integer; number of initial samples to discard from each
#'   chain.
#' @param max_iter Numeric; maximum iteration to retain (samples beyond this
#'   are dropped).
#' @return A `multiPhylo` object containing the pooled post-burn-in trees.
#' @export
collect_chains = function(edge_list_all, phy_init, burnin = 0, max_iter = Inf) {

    if (max_iter < burnin) {
        stop('Max iter needs to be greater than burnin')
    }

    # drop empty chains
    edge_list_all = edge_list_all[edge_list_all %>% sapply(length) > 0]

    # reconstruct trees per chain and apply burnin/truncation
    mcmc_trees = edge_list_all %>% lapply(function(elist){
            elist = elist[(burnin+1):min(length(elist), max_iter)]
            lapply(elist, function(edges){
                tree = attach_edges(phy_init, edges)
                tree$edge.length = NULL
                tree$node.label = NULL
                tree$nodes = NULL
                tree$edge = TreeTools::RenumberTree(tree$edge[,1], tree$edge[,2])
                return(tree)
            })
        }) %>% unlist(recursive = F)

    class(mcmc_trees) = 'multiPhylo'

    return(mcmc_trees)
    
}

#' Collect raw edge matrices from MCMC chains
#'
#' Pools edge matrices from all chains after applying burn-in removal and
#' iteration truncation. Unlike [collect_chains()], does not reconstruct
#' full `phylo` objects.
#'
#' @param edge_list_all List of chains, each a list of edge matrices.
#' @param burnin Integer; number of initial samples to discard.
#' @param max_iter Numeric; maximum iteration to retain.
#' @return A flat list of 2-column edge matrices.
#' @keywords internal
#' @noRd
collect_edges = function(edge_list_all, burnin = 0, max_iter = Inf) {

    if (max_iter < burnin) {
        stop('Max iter needs to be greater than burnin')
    }

    mcmc_edges = edge_list_all %>% lapply(function(elist){
            elist = elist[(burnin+1):min(length(elist), max_iter)]
        }) %>% unlist(recursive = F)

    return(mcmc_edges)
}

#' Add clade frequencies to a phylogenetic tree
#'
#' Computes the frequency (support) of each clade in a reference phylogeny (`phy`) across a list of phylogenetic trees (`edge_list`),
#' and adds these frequencies as node labels to the reference tree.
#'
#' @param phy A reference phylogeny of class \code{phylo}.
#' @param edge_list A list of edge matrices (or phylogenetic trees) to compare against the reference tree.
#' @param rooted Logical; whether to treat the trees as rooted. Default is \code{TRUE}.
#' @param ncores Integer; number of cores to use for parallel computation. Default is \code{1} (no parallelization).
#' @return A phylo object with clade frequencies added as node labels.
#' @keywords internal
add_clade_freq = function(phy, edge_list, rooted = TRUE, ncores = 1) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)
    phy = reorder_phylo(phy) # prop_clades_par requires phylo in postorder
    freqs = prop_clades_par(phy$edge, edge_list, rooted = rooted, normalize = TRUE)
    phy$node.label = freqs
    return(phy)
}

#' Add clade frequencies to a tree by streaming an on-disk trace
#'
#' Equivalent to `add_clade_freq(phy, collect_edges(safe_read_chain(path), ...))`,
#' but never materializes the whole trace in R at once. `safe_read_chain()`
#' reads every block file (and the legacy file, if any) fully into memory
#' before burn-in is applied -- for a long block-wise trace (one block per
#' MCMC batch, `run_tree_mcmc_batch()`'s default `checkpoint_every = 1`) the
#' decompressed in-memory size can be an order of magnitude larger than the
#' on-disk (qs2/zstd-compressed) size, since raw integer edge matrices
#' compress heavily but every element still has to become a live R object to
#' be scanned. `prop_clades_par()` only needs to *scan* trees to accumulate
#' bipartition counts (`normalize = FALSE`), which is associative across any
#' partition of the tree list -- so this reads one block at a time, trims
#' burnin/max_iter using a running per-chain element count (reproducing
#' exactly the slice `collect_edges()` would take from the fully-concatenated
#' list), accumulates raw counts, and discards the block before reading the
#' next. Peak memory is bounded by one block (`checkpoint_every` batches x
#' `nchains`), not by the whole trace.
#'
#' @param phy A reference phylogeny of class `phylo` (the base tree to annotate).
#' @param path Trace file path, as passed to `safe_read_chain()`/`run_tree_mcmc_batch()`.
#' @param burnin Integer; number of initial samples to discard per chain.
#' @param max_iter Numeric; maximum per-chain sample index to retain.
#' @param rooted Logical; whether to treat the trees as rooted. Default `TRUE`.
#' @param ncores Integer; threads for the clade-counting scan (`RcppParallel`).
#' @param ncores_qs Integer; threads for `qs2` deserialization of each block.
#' @return A phylo object with clade frequencies added as node labels.
#' @keywords internal
#' @noRd
streaming_clade_freq = function(phy, path, burnin = 0, max_iter = Inf, rooted = TRUE, ncores = 1, ncores_qs = 1) {
    if (max_iter < burnin) {
        stop('Max iter needs to be greater than burnin')
    }
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)
    phy = reorder_phylo(phy) # prop_clades_par requires phylo in postorder

    legacy_exists = file.exists(path) && !dir.exists(path)
    block_files = list_trace_blocks(path)
    if (!legacy_exists && length(block_files) == 0) {
        stop('No MCMC trace found at ', path, ' (or its .blocks/ directory)')
    }

    counts = NULL
    total_n = 0
    seen = list() # per-chain count of elements consumed so far, across chunks in order

    consume_chunk = function(chunk) {
        keep = list()
        for (nm in names(chunk)) {
            elist = chunk[[nm]]
            L = length(elist)
            if (L == 0) next
            prev_seen = if (is.null(seen[[nm]])) 0L else seen[[nm]]
            seen[[nm]] <<- prev_seen + L
            # local (within-chunk) bounds equivalent to the global
            # [burnin+1, max_iter] slice collect_edges() would take from the
            # fully-concatenated per-chain list
            lo = max(1L, burnin - prev_seen + 1L)
            hi = min(L, max_iter - prev_seen)
            if (lo > hi) next
            keep = c(keep, elist[lo:hi])
        }
        if (length(keep) == 0) return(invisible())
        raw = prop_clades_par(phy$edge, keep, rooted = rooted, normalize = FALSE)
        counts <<- if (is.null(counts)) raw else counts + raw
        total_n <<- total_n + length(keep)
    }

    if (legacy_exists) {
        legacy = read_qs2_file(path, ncores = ncores_qs)
        if (!is.null(legacy)) consume_chunk(legacy)
        rm(legacy)
    }
    for (bf in block_files) {
        blk = read_qs2_file(bf, ncores = ncores_qs)
        if (!is.null(blk)) consume_chunk(blk)
        rm(blk)
    }

    if (is.null(counts) || total_n == 0) {
        stop('No post-burnin MCMC samples available to compute clade frequencies')
    }
    phy$node.label = counts / total_n
    return(phy)
}
