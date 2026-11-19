# Branch-length (node-height) inference on a fixed topology.
#
# Each edge gets its own Wright-Fisher transition over t_e = h(parent) - h(child) generations, interpolated between
# the matrices of a precomputed grid (as in the original branch-length work, interpolate_matrices()). Mutation
# scales with edge length: the package's per-edge eps (fitted for an edge of `ngen_ref` generations) is converted
# to a per-generation rate, so an edge of t generations uses eps_t = 1 - (1 - eps_gen)^t and a zero-length edge
# creates no new variant. The C++ sampler (src/brl.cpp) updates only the edges at the moved node and the messages
# on its path to the root.

#' Transition-matrix grid for branch-length inference
#'
#' @param k Number of interior VAF bins (matrices are `(k + 2) x (k + 2)`).
#' @param npop mtDNA copies per cell (`N` of the WF model).
#' @param eps Mutation (boundary) probability per edge of `ngen_ref` generations, as fitted by the package.
#' @param ngen_ref Generations per edge that `eps` refers to (the fitted `ngen`).
#' @param grid Increasing generation values (integers, starting at 0).
#' @return A list with `grid` and `A` (row-major probability vectors, one per grid point), `eps_gen`.
#' @export
brl_transition_grid <- function(k, npop, eps, ngen_ref, grid) {
    grid <- sort(unique(as.integer(round(grid))))
    if (grid[1] != 0) grid <- c(0L, grid)
    eps_gen <- 1 - (1 - eps)^(1 / ngen_ref)
    A <- lapply(grid, function(t) {
        At <- get_transition_mat_wf_hmm(k = k, eps = 1 - (1 - eps_gen)^t, N = npop, ngen = t, safe = TRUE)
        as.vector(t(At))
    })
    list(grid = as.numeric(grid), A = A, eps_gen = eps_gen, k = k, npop = npop)
}

#' Node heights of a phylo object (tips at 0)
#'
#' @param phy A rooted `phylo` with `edge.length` (made ultrametric by taking each tip's depth to the deepest tip).
#' @return Numeric vector of heights indexed by node number.
#' @export
brl_heights_from_phylo <- function(phy) {
    d <- ape::node.depth.edgelength(phy)
    max(d[seq_along(phy$tip.label)]) - d
}

#' Grafen-style starting heights
#'
#' Node height = `htotal * (descendant tips - 1) / (all tips - 1)` (root at `htotal`, tips at 0), so every parent is
#' above its children.
#' @param phy A rooted `phylo`.
#' @param htotal Root height in generations.
#' @export
brl_init_heights <- function(phy, htotal) {
    nt <- length(phy$tip.label)
    E <- ape::reorder.phylo(phy, "postorder")$edge
    ndesc <- c(rep(1, nt), rep(0, phy$Nnode))
    for (i in seq_len(nrow(E))) ndesc[E[i, 1]] <- ndesc[E[i, 1]] + ndesc[E[i, 2]]
    h <- htotal * (ndesc - 1) / (nt - 1)
    h[nt + 1] <- htotal
    h
}

brl_setup <- function(phy, logliks) {
    phy$node.label <- NULL
    E <- reorder_phylo(phy)$edge
    list(E = as.integer(c(E[, 1], E[, 2])), logP = convert_logliks_to_logP_list(logliks, phy), phy = phy)
}

#' Log-likelihood of a tree with branch lengths
#'
#' @param phy Rooted `phylo` with node heights given by `heights` (or by its `edge.length` if `heights` is NULL).
#' @param logliks Leaf log-likelihoods (`md$leaf_likelihoods`).
#' @param grid_obj Output of [brl_transition_grid()].
#' @param heights Optional node heights (tips 0) indexed by node number.
#' @param ncores Threads.
#' @return Total log-likelihood; per-locus values in attribute `per_locus`.
#' @export
score_tree_brl <- function(phy, logliks, grid_obj, heights = NULL, ncores = 1) {
    RcppParallel::setThreadOptions(numThreads = ncores)
    s <- brl_setup(phy, logliks)
    if (is.null(heights)) heights <- brl_heights_from_phylo(phy)
    ll <- score_tree_brl_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, as.numeric(heights))
    structure(sum(ll), per_locus = setNames(ll, names(s$logP)))
}

#' MCMC over node heights on a fixed topology
#'
#' Random-walk Metropolis-Hastings on each internal node's position between its highest child and its parent
#' (logit scale), root height fixed at `htotal`, uniform prior on valid heights. Only the moved node's edges and its
#' path to the root are rescored per proposal.
#'
#' @param phy Rooted `phylo` (topology used; polytomies allowed).
#' @param logliks Leaf log-likelihoods (`md$leaf_likelihoods`), named by variant, columns = cells.
#' @param k,npop,eps,ngen_ref Model parameters (`md$model_params`: k, npop, eps, ngen).
#' @param htotal Root height in generations (fixed).
#' @param n_iter Iterations (one node move each).
#' @param grid Generation grid for the transition matrices (default: 0 to `htotal`, about 100 points).
#' @param heights_init Starting heights (default: Grafen-style).
#' @param step_sd Random-walk sd on the logit scale.
#' @param thin Keep every `thin`-th state.
#' @param burnin Fraction of kept states dropped before averaging.
#' @param seed RNG seed.
#' @param ncores Threads (parallel over variants).
#' @param refresh Full likelihood recomputation every `refresh` iterations.
#' @return A list: `trace` (kept heights, one column per internal non-root node), `logL`, `acceptance`,
#'   `heights_mean` (posterior mean, all nodes), `phy_mean` (`phy` with branch lengths from the posterior mean),
#'   `logL_check` (incremental vs full recomputation at the end), `grid_obj`.
#' @export
run_brl_mcmc <- function(phy, logliks, k, npop, eps, ngen_ref, htotal, n_iter = 10000, grid = NULL,
                         heights_init = NULL, step_sd = 1, thin = 10, burnin = 0.25, seed = 1, ncores = 1,
                         refresh = 1000) {
    RcppParallel::setThreadOptions(numThreads = ncores)
    if (is.null(grid)) grid <- unique(round(seq(0, ceiling(htotal), length.out = 101)))
    grid_obj <- brl_transition_grid(k = k, npop = npop, eps = eps, ngen_ref = ngen_ref, grid = grid)
    s <- brl_setup(phy, logliks)
    if (is.null(heights_init)) heights_init <- brl_init_heights(s$phy, htotal)
    res <- tree_brl_mcmc_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, as.numeric(heights_init), n_iter = n_iter,
                             step_sd = step_sd, thin = thin, seed = seed, refresh = refresh)
    keep <- seq_len(nrow(res$heights)); keep <- keep[keep > floor(burnin * length(keep))]
    hm <- as.numeric(heights_init)
    hm[res$nodes] <- colMeans(res$heights[keep, , drop = FALSE])
    phy_mean <- s$phy; phy_mean$edge.length <- hm[phy_mean$edge[, 1]] - hm[phy_mean$edge[, 2]]
    colnames(res$heights) <- res$nodes
    list(trace = res$heights, logL = res$logL, accepted = res$accepted, acceptance = res$acceptance,
         nodes = res$nodes, heights_mean = hm, heights_final = res$heights_final, phy_mean = phy_mean,
         logL_check = c(incremental = res$logL_incremental, full = res$logL_full), grid_obj = grid_obj,
         htotal = htotal)
}

#' Branch-length MCMC from a fitted MitoDrift object
#'
#' @param md A fitted `MitoDrift` object.
#' @param phy Topology (default `md$tree_ml`, else `md$tree_init`).
#' @param htotal Root height in generations (default: fitted `ngen` times the mean number of edges from the root to
#'   a tip, i.e. the total drift the fixed-length model implies).
#' @param ... Passed to [run_brl_mcmc()].
#' @export
md_brl_mcmc <- function(md, phy = NULL, htotal = NULL, ...) {
    if (is.null(phy)) phy <- if (!is.null(md$tree_ml)) md$tree_ml else md$tree_init
    if (inherits(phy, "multiPhylo")) phy <- phy[[length(phy)]]
    p <- md$model_params
    if (is.null(htotal)) {
        nt <- length(phy$tip.label)
        depth_edges <- ape::node.depth.edgelength(ape::compute.brlen(phy, 1))[seq_len(nt)]
        htotal <- unname(p["ngen"]) * mean(depth_edges)
    }
    run_brl_mcmc(phy, md$leaf_likelihoods, k = unname(p["k"]), npop = unname(p["npop"]), eps = unname(p["eps"]),
                 ngen_ref = unname(p["ngen"]), htotal = htotal, ...)
}
