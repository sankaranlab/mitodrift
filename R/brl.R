# Branch-length (node-height) inference on a fixed topology.
#
# Each edge gets its own Wright-Fisher transition over t_e = h(parent) - h(child) generations, interpolated between
# the matrices of a precomputed grid (as in the original branch-length work, interpolate_matrices()). Mutation
# scales with edge length: the package's per-edge eps (fitted for an edge of `ngen_ref` generations) is converted
# to a per-generation rate, so an edge of t generations uses eps_t = 1 - (1 - eps_gen)^t and a zero-length edge
# creates no new variant. The C++ sampler (src/brl.cpp) updates only the edges at the moved node and the messages
# on its path to the root. The root height is sampled by a global scaling move (all heights x c; 1/H prior).

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

# Edge vector and per-locus logP for the C++ code. Uses ape's postorder (root = Ntip + 1, root edges last), not
# reorder_phylo()/reorderRcpp(), which assume a binary tree (n_tips = n_edges / 2 + 1); polytomies are allowed here.
brl_setup <- function(phy, logliks) {
    phy$node.label <- NULL
    E <- ape::reorder.phylo(phy, "postorder")$edge
    nt <- length(phy$tip.label); nn <- phy$Nnode; root <- nt + 1L
    logP <- lapply(logliks, function(lm) {
        k <- nrow(lm); P <- matrix(log(1 / k), k, nt + nn)
        P[, seq_len(nt)] <- lm[, phy$tip.label]
        P[, root] <- log(c(1, rep(0, k - 1)))
        as.vector(t(P))
    })
    list(E = as.integer(c(E[, 1], E[, 2])), logP = logP, phy = phy)
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
#' @param htotal Starting root height in generations (fixed if `sample_root = FALSE`).
#' @param n_iter Iterations (one node move each).
#' @param grid Generation grid for the transition matrices (default: 0 to `htotal`, about 100 points).
#' @param heights_init Starting heights (default: Grafen-style).
#' @param step_sd Random-walk sd on the logit scale.
#' @param thin Keep every `thin`-th state.
#' @param burnin Fraction of kept states dropped before averaging.
#' @param seed RNG seed.
#' @param ncores Threads (parallel over variants).
#' @param refresh Full likelihood recomputation every `refresh` iterations.
#' @param sample_root Sample the root height (global scaling move, 1/H prior); default TRUE.
#' @param p_scale,scale_sd Probability of a scaling move per iteration and sd of its log scale factor.
#' @param hmax_factor Grid (and root height) upper bound as a multiple of `htotal`.
#' @return A list: `trace` (kept heights in generations, one column per internal non-root node), `trace_rel` (the same
#'   divided by each sample's root height), `root_height`, `logL`, `acceptance`, `heights_rel_mean` (posterior mean
#'   relative heights, all nodes), `root_height_mean`, `heights_mean` (relative mean x root mean), `phy_mean`,
#'   `logL_check` (incremental vs full recomputation at the end), `grid_obj`, `keep` (post-burn-in rows).
#' @export
run_brl_mcmc <- function(phy, logliks, k, npop, eps, ngen_ref, htotal, n_iter = 10000, grid = NULL,
                         heights_init = NULL, step_sd = 1, thin = 10, burnin = 0.25, seed = 1, ncores = 1,
                         refresh = 1000, sample_root = TRUE, p_scale = 0.02, scale_sd = 0.15, hmax_factor = 50) {
    RcppParallel::setThreadOptions(numThreads = ncores)
    if (is.null(grid)) {
        hmax <- if (sample_root) ceiling(htotal * hmax_factor) else ceiling(htotal)
        # integer generations, dense near 0 (short edges) and log-spaced above, so any root height up to hmax is covered
        grid <- sort(unique(c(0:min(30, hmax), round(exp(seq(log(min(31, hmax)), log(hmax), length.out = 90))))))
    }
    grid_obj <- brl_transition_grid(k = k, npop = npop, eps = eps, ngen_ref = ngen_ref, grid = grid)
    s <- brl_setup(phy, logliks)
    if (is.null(heights_init)) heights_init <- brl_init_heights(s$phy, htotal)
    res <- tree_brl_mcmc_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, as.numeric(heights_init), n_iter = n_iter,
                             step_sd = step_sd, thin = thin, seed = seed, refresh = refresh,
                             p_scale = if (sample_root) p_scale else 0, scale_sd = scale_sd)
    keep <- seq_len(nrow(res$heights)); keep <- keep[keep > floor(burnin * length(keep))]
    colnames(res$heights) <- res$nodes
    rel <- res$heights / res$root_height                     # heights relative to each sample's root height
    root <- nt_root <- length(s$phy$tip.label) + 1
    H_mean <- mean(res$root_height[keep])
    hrel <- c(rep(0, length(s$phy$tip.label)), rep(NA, s$phy$Nnode)); hrel[root] <- 1
    hrel[res$nodes] <- colMeans(rel[keep, , drop = FALSE])
    hm <- hrel * H_mean                                      # posterior mean shape at the posterior mean root height
    phy_mean <- s$phy; phy_mean$edge.length <- hm[phy_mean$edge[, 1]] - hm[phy_mean$edge[, 2]]
    list(trace = res$heights, trace_rel = rel, root_height = res$root_height, logL = res$logL, accepted = res$accepted,
         acceptance = res$acceptance, scale_acceptance = res$scale_acceptance, nodes = res$nodes,
         heights_mean = hm, heights_rel_mean = hrel, root_height_mean = H_mean, heights_final = res$heights_final,
         phy_mean = phy_mean, logL_check = c(incremental = res$logL_incremental, full = res$logL_full),
         grid_obj = grid_obj, htotal_init = htotal, keep = keep)
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


#' Maximum-likelihood node heights on a fixed topology
#'
#' `method = "joint"` (default): all node positions and the root height are optimised together by L-BFGS-B, in the
#' top-down parameterisation of the original `x_to_length()`: a node's depth from the root is
#' `d_v = d_parent + (H - d_parent) * x_v`, `x_v = sigmoid(p_v)`, so moving a node carries its subtree and every
#' parameter vector is a valid ultrametric tree. The likelihood is a full C++ pass; the gradient is exact, from one
#' inside-outside pass ([brl_grad_exact_cpp], parallel over loci), chained analytically into `p` and log H. `method = "coordinate"`: per-node golden-section sweeps (as the original
#' `estimate_brl_iter()`), each node moving between its highest child and its parent.
#'
#' @inheritParams run_brl_mcmc
#' @param method "joint" or "coordinate".
#' @param optimize_root Also optimise the root height.
#' @param max_sweeps,tol_ll Stopping rule for "coordinate"; `maxit` for "joint".
#' @return A list: `heights` (generations, indexed by node), `heights_rel` (root = 1), `root_height`, `logL`, `phy`
#'   (with ML branch lengths), `convergence` (optim code for "joint"), `n_eval`.
#' @export
optimize_brl <- function(phy, logliks, k, npop, eps, ngen_ref, htotal, grid = NULL, heights_init = NULL,
                         method = c("joint", "coordinate"), optimize_root = TRUE, max_sweeps = 50, tol_ll = 0.01,
                         maxit = 500, hmax_factor = 50, ncores = 1) {
    method <- match.arg(method)
    RcppParallel::setThreadOptions(numThreads = ncores)
    if (is.null(grid)) {
        hmax <- if (optimize_root) ceiling(htotal * hmax_factor) else ceiling(htotal)
        grid <- sort(unique(c(0:min(30, hmax), round(exp(seq(log(min(31, hmax)), log(hmax), length.out = 90))))))
    }
    grid_obj <- brl_transition_grid(k = k, npop = npop, eps = eps, ngen_ref = ngen_ref, grid = grid)
    s <- brl_setup(phy, logliks)
    if (is.null(heights_init)) heights_init <- brl_init_heights(s$phy, htotal)
    nt <- length(s$phy$tip.label); root <- nt + 1L
    if (method == "coordinate") {
        r <- tree_brl_optimize_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, as.numeric(heights_init), max_sweeps = max_sweeps,
                                   tol_ll = tol_ll, optimize_root = optimize_root, max_factor = hmax_factor)
        h <- r$heights; conv <- NA_integer_; n_eval <- r$sweeps; ll <- r$logL
    } else {
        Ec <- ape::reorder.phylo(s$phy, "cladewise")$edge          # parents before children
        par_of <- integer(nt + s$phy$Nnode); par_of[Ec[, 2]] <- Ec[, 1]
        inner <- Ec[Ec[, 2] > nt, 2]                                # internal non-root nodes, top-down
        inner_po <- rev(inner)                                      # bottom-up
        kids <- split(Ec[, 2], Ec[, 1])
        to_heights <- function(pv, H) {
            d <- numeric(nt + s$phy$Nnode); d[seq_len(nt)] <- H; x <- 1 / (1 + exp(-pv))
            for (i in seq_along(inner)) { v <- inner[i]; dp <- d[par_of[v]]; d[v] <- dp + (H - dp) * x[i] }
            list(h = H - d, d = d, x = x)
        }
        H0 <- heights_init[root]; d0 <- H0 - heights_init
        x0 <- (d0[inner] - d0[par_of[inner]]) / (H0 - d0[par_of[inner]]); x0 <- pmin(pmax(x0, 1e-4), 1 - 1e-4)
        par0 <- c(log(x0 / (1 - x0)), if (optimize_root) log(H0))
        np <- length(inner); n_eval <- 0L
        getH <- function(par) if (optimize_root) exp(par[np + 1]) else H0
        fn <- function(par) { n_eval <<- n_eval + 1L; z <- to_heights(par[seq_len(np)], getH(par))
            -sum(score_tree_brl_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, z$h)) }
        gr <- function(par) {
            H <- getH(par); z <- to_heights(par[seq_len(np)], H)
            g <- brl_grad_exact_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, z$h)          # exact, one inside-outside pass
            gh <- numeric(nt + s$phy$Nnode); gh[as.integer(names(g$grad))] <- g$grad
            gd <- -gh; gd[root] <- 0                                                         # dlogL/d(depth) = -dlogL/dh
            adj <- gd; xi <- setNames(z$x, inner)
            for (v in inner_po) for (c in kids[[as.character(v)]]) if (c > nt) adj[v] <- adj[v] + adj[c] * (1 - xi[[as.character(c)]])
            gp <- vapply(seq_along(inner), function(i) { v <- inner[i]; adj[v] * (H - z$d[par_of[v]]) * z$x[i] * (1 - z$x[i]) }, 0)
            out <- -gp
            if (optimize_root) out <- c(out, -sum(gh * z$h))   # scaling every height by c: dlogL/dlog c = sum_v h_v dlogL/dh_v
            out
        }
        lo <- c(rep(-15, np), if (optimize_root) log(max(H0 / hmax_factor, 1e-3)))
        hi <- c(rep(15, np), if (optimize_root) log(max(grid_obj$grid)))
        o <- stats::optim(par0, fn, gr, method = "L-BFGS-B", lower = lo, upper = hi, control = list(maxit = maxit))
        z <- to_heights(o$par[seq_len(np)], getH(o$par)); h <- z$h; conv <- o$convergence; ll <- -o$value
    }
    phy_ml <- s$phy; phy_ml$edge.length <- h[phy_ml$edge[, 1]] - h[phy_ml$edge[, 2]]
    list(heights = h, heights_rel = h / h[root], root_height = h[root], logL = ll, phy = phy_ml, convergence = conv,
         n_eval = n_eval, method = method, grid_obj = grid_obj)
}

#' Hamiltonian Monte Carlo over node heights on a fixed topology
#'
#' Samples the same posterior as [run_brl_mcmc()] (uniform prior on valid node heights given the root height H, prior
#' 1/H on H) in the unconstrained top-down parameterisation of `x_to_length()`: a node's relative depth is
#' `D_v = D_parent + (1 - D_parent) x_v`, `x_v = sigmoid(p_v)`, heights `h = H (1 - D)`, plus `log H`. In these
#' coordinates the prior with its Jacobian is `prod_v (1 - D_parent(v)) x_v (1 - x_v)`, whose log-gradient is
#' `(1 - 2 x_a) - n_a x_a` for node a with n_a internal descendants. The likelihood gradient is exact (one
#' inside-outside pass per locus, [brl_grad_exact_cpp()], parallel over loci), chained to `p` and `log H`.
#' Every iteration moves all nodes (leapfrog integration, `n_leapfrog` steps); warm-up tunes the step size by dual
#' averaging (target acceptance `delta`) and a diagonal mass matrix from the middle of warm-up.
#'
#' @inheritParams run_brl_mcmc
#' @param n_warmup,n_sample Warm-up and kept iterations.
#' @param n_leapfrog Leapfrog steps per iteration (jittered +-50%).
#' @param delta Target acceptance rate.
#' @param heights_init Starting heights (default: Grafen-style at `htotal`).
#' @return A list: `trace_rel` (kept relative heights, one column per internal non-root node), `root_height`, `logL`,
#'   `accept` (per kept iteration), `step_size`, `nodes`, `seconds`.
#' @export
run_brl_hmc <- function(phy, logliks, k, npop, eps, ngen_ref, htotal, n_warmup = 500, n_sample = 1000, n_leapfrog = 20,
                        delta = 0.8, grid = NULL, heights_init = NULL, hmax_factor = 50, seed = 1, ncores = 1) {
    set.seed(seed); t0 <- Sys.time()
    RcppParallel::setThreadOptions(numThreads = ncores)
    if (is.null(grid)) {
        hmax <- ceiling(htotal * hmax_factor)
        grid <- sort(unique(c(0:min(30, hmax), round(exp(seq(log(min(31, hmax)), log(hmax), length.out = 90))))))
    }
    grid_obj <- brl_transition_grid(k = k, npop = npop, eps = eps, ngen_ref = ngen_ref, grid = grid)
    s <- brl_setup(phy, logliks)
    nt <- length(s$phy$tip.label); nn <- s$phy$Nnode; root <- nt + 1L
    if (is.null(heights_init)) heights_init <- brl_init_heights(s$phy, htotal)
    Ec <- ape::reorder.phylo(s$phy, "cladewise")$edge
    par_of <- integer(nt + nn); par_of[Ec[, 2]] <- Ec[, 1]
    inner <- Ec[Ec[, 2] > nt, 2]; np <- length(inner); pos <- integer(nt + nn); pos[inner] <- seq_len(np)
    kids <- split(Ec[, 2], Ec[, 1]); inner_po <- rev(inner)
    ndesc <- integer(nt + nn); for (v in inner_po) for (c in kids[[as.character(v)]]) if (c > nt) ndesc[v] <- ndesc[v] + 1L + ndesc[c]
    hmax_val <- max(grid_obj$grid)
    geom <- function(th) {
        x <- 1 / (1 + exp(-th[seq_len(np)])); H <- exp(th[np + 1])
        D <- numeric(nt + nn); D[seq_len(nt)] <- 1
        for (i in seq_len(np)) { v <- inner[i]; Dp <- D[par_of[v]]; D[v] <- Dp + (1 - Dp) * x[i] }
        list(x = x, H = H, D = D, h = H * (1 - D))
    }
    logpost_grad <- function(th) {
        g <- geom(th)
        if (g$H > hmax_val) return(list(lp = -Inf, grad = rep(0, np + 1)))
        r <- brl_grad_exact_cpp(s$E, s$logP, grid_obj$A, grid_obj$grid, g$h)
        gh <- numeric(nt + nn); gh[as.integer(names(r$grad))] <- r$grad
        # likelihood: dlogL/dD_v = -H dlogL/dh_v (non-root); adjoint down the tree; then to p and log H
        adj <- -g$H * gh; adj[root] <- 0
        for (v in inner_po) for (c in kids[[as.character(v)]]) if (c > nt) adj[v] <- adj[v] + adj[c] * (1 - g$x[pos[c]])
        Dp <- g$D[par_of[inner]]
        gp <- adj[inner] * (1 - Dp) * g$x * (1 - g$x)
        glogH <- sum(gh * g$h)
        # prior with Jacobian
        lp_prior <- sum(log1p(-Dp)) + sum(log(g$x) + log1p(-g$x))
        gp <- gp + (1 - 2 * g$x) - ndesc[inner] * g$x
        list(lp = r$logL + lp_prior, grad = c(gp, glogH), logL = r$logL, H = g$H, h = g$h)
    }
    # start
    H0 <- heights_init[root]; D0 <- 1 - heights_init / H0
    x0 <- (D0[inner] - D0[par_of[inner]]) / pmax(1 - D0[par_of[inner]], 1e-12); x0 <- pmin(pmax(x0, 1e-4), 1 - 1e-4)
    th <- c(log(x0 / (1 - x0)), log(H0)); cur <- logpost_grad(th)
    minv <- rep(1, np + 1)                                   # inverse mass (variances)
    # step size: find a reasonable start, then dual averaging during warm-up
    eps_s <- 0.01; mu <- log(10 * eps_s); Hbar <- 0; logeps_bar <- 0; gamma <- 0.05; t0da <- 10; kap <- 0.75
    n_tot <- n_warmup + n_sample; keep_rows <- (n_warmup + 1):n_tot
    trace <- matrix(NA_real_, n_sample, np); root_tr <- numeric(n_sample); ll_tr <- numeric(n_sample); acc_tr <- numeric(n_sample)
    warm_store <- matrix(NA_real_, 0, np + 1); win <- c(floor(0.3 * n_warmup), floor(0.85 * n_warmup))
    for (it in seq_len(n_tot)) {
        L <- max(1L, round(n_leapfrog * runif(1, 0.5, 1.5)))
        r0 <- rnorm(np + 1) / sqrt(minv); th1 <- th; r <- r0 + 0.5 * eps_s * cur$grad; nw <- cur; ok <- TRUE
        for (l in seq_len(L)) {
            th1 <- th1 + eps_s * minv * r
            nw <- logpost_grad(th1); if (!is.finite(nw$lp)) { ok <- FALSE; break }
            if (l < L) r <- r + eps_s * nw$grad
        }
        if (ok) r <- r + 0.5 * eps_s * nw$grad
        logacc <- if (ok) (nw$lp - 0.5 * sum(r^2 * minv)) - (cur$lp - 0.5 * sum(r0^2 * minv)) else -Inf
        a <- min(1, exp(logacc)); if (!is.finite(a)) a <- 0
        if (runif(1) < a) { th <- th1; cur <- nw }
        if (it <= n_warmup) {
            m <- it; Hbar <- (1 - 1 / (m + t0da)) * Hbar + (delta - a) / (m + t0da)
            logeps <- mu - sqrt(m) / gamma * Hbar; eta <- m^(-kap); logeps_bar <- eta * logeps + (1 - eta) * logeps_bar
            eps_s <- exp(logeps)
            if (it > win[1] && it <= win[2]) warm_store <- rbind(warm_store, th)
            if (it == win[2] && nrow(warm_store) > 20) {        # diagonal mass matrix, then restart step-size adaptation
                minv <- pmax(apply(warm_store, 2, var), 1e-6); mu <- log(10 * eps_s); Hbar <- 0; logeps_bar <- 0
            }
            if (it == n_warmup) eps_s <- exp(logeps_bar)
        } else {
            j <- it - n_warmup; g <- geom(th)
            trace[j, ] <- (g$h / g$H)[inner]; root_tr[j] <- g$H; ll_tr[j] <- cur$logL; acc_tr[j] <- a
        }
    }
    colnames(trace) <- inner
    list(trace_rel = trace, root_height = root_tr, logL = ll_tr, accept = acc_tr, step_size = eps_s, nodes = inner,
         minv = minv, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
