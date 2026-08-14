make_mc3_test_model <- function(n_tips = NULL) {
    mut_dat <- utils::read.csv(system.file("extdata", "small_test_mut_dat.csv", package = "mitodrift"))
    if (!is.null(n_tips)) {
        cells <- sort(unique(mut_dat$cell))
        stopifnot(n_tips >= 3L, n_tips <= length(cells))
        mut_dat <- mut_dat[mut_dat$cell %in% cells[seq_len(n_tips)], , drop = FALSE]
    }
    md <- MitoDrift$new(
        mut_dat = mut_dat,
        model_params = c(eps = 0.001, err = 0, npop = 600, ngen = 100, k = 20),
        ncores = 1
    )
    md$make_model(ncores = 1)
    md
}

rooted_topology_key <- function(edge, n_tips, reorder = FALSE) {
    if (reorder) edge <- mitodrift:::reorderRcpp(edge)
    edge <- matrix(edge, ncol = 2L)
    stopifnot(n_tips <= 30L)

    descendant_bits <- integer(max(edge))
    descendant_bits[seq_len(n_tips)] <- bitwShiftL(1L, seq_len(n_tips) - 1L)
    for (i in seq_len(nrow(edge))) {
        parent <- edge[i, 1L]
        child <- edge[i, 2L]
        descendant_bits[parent] <- bitwOr(descendant_bits[parent], descendant_bits[child])
    }

    root <- edge[nrow(edge), 1L]
    clades <- descendant_bits[setdiff(unique(edge[, 1L]), root)]
    paste(sort(clades), collapse = "-")
}

exact_rooted_topology_posterior <- function(md) {
    n_tips <- length(md$tree_init$tip.label)
    trees <- phangorn::allTrees(
        n_tips, rooted = TRUE, tip.label = md$tree_init$tip.label)
    keys <- vapply(
        trees, function(tree) rooted_topology_key(tree$edge, n_tips, reorder = TRUE),
        character(1L))
    stopifnot(!anyDuplicated(keys))
    log_likelihood <- vapply(trees, function(tree) {
        mitodrift:::score_tree_bp_wrapper2(tree$edge, md$logP, md$logA)
    }, numeric(1L))
    probability <- exp(log_likelihood - max(log_likelihood))
    probability <- probability / sum(probability)
    names(probability) <- keys

    list(
        trees = trees,
        keys = keys,
        log_likelihood = log_likelihood,
        probability = probability
    )
}

estimate_topology_posterior <- function(traces, keys, n_tips, burnin) {
    chain_estimates <- lapply(traces, function(chain) {
        stopifnot(length(chain) > burnin + 1L)
        sampled_keys <- vapply(
            chain[(burnin + 2L):length(chain)],
            rooted_topology_key,
            character(1L),
            n_tips = n_tips
        )
        stopifnot(all(sampled_keys %in% keys))
        counts <- table(factor(sampled_keys, levels = keys))
        as.numeric(counts) / length(sampled_keys)
    })
    rowMeans(do.call(cbind, chain_estimates))
}
