test_that("one-temperature DEO MC3 exactly matches ordinary MCMC", {
    md <- make_mc3_test_model()
    niter <- 100L
    seed <- 42L

    ordinary <- mitodrift:::tree_mcmc_cpp_cached_threadsafe(
        md$tree_init$edge, md$logP, md$logA, niter, seed, TRUE)
    deo <- mitodrift:::tree_mc3_parallel_seeded_deo(
        list(list(md$tree_init$edge)), md$logP, md$logA,
        niter, seed, 1, 10L)

    expect_identical(deo$traces[[1]], ordinary)
    expect_equal(dim(deo$swap_attempts), c(1L, 0L))
    expect_equal(dim(deo$swap_accepts), c(1L, 0L))
})

test_that("DEO alternates every disjoint adjacent pair with exact attempt counts", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.05, 1.1, 1.15, 1.2)
    starts <- rep(
        list(rep(list(md$tree_init$edge), length(temperatures))),
        2L)

    deo <- mitodrift:::tree_mc3_parallel_seeded_deo(
        starts, md$logP, md$logA, c(25L, 16L), c(7L, 8L),
        temperatures, 5L)

    # Five and three complete swap barriers, respectively. Even-indexed
    # pairs are attempted at barriers 1,3,5; odd-indexed pairs at 2,4.
    expect_identical(
        deo$swap_attempts,
        matrix(c(3L, 2L, 3L, 2L,
                 2L, 1L, 2L, 1L), nrow = 2L, byrow = TRUE))
    expect_true(all(deo$swap_accepts >= 0L))
    expect_true(all(deo$swap_accepts <= deo$swap_attempts))
    expect_equal(lengths(deo$traces), c(26L, 17L))
    expect_equal(lengths(deo$final_states), c(5L, 5L))

    scores <- unlist(lapply(deo$final_states, function(states) {
        vapply(states, function(edges) {
            mitodrift:::score_tree_bp_wrapper2(edges, md$logP, md$logA)
        }, numeric(1))
    }))
    expect_true(all(is.finite(scores)))
})

test_that("DEO records the cold state after the final swap barrier", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.0001, 1.0002, 1.0003)
    niter <- 200L

    deo <- mitodrift:::tree_mc3_parallel_seeded_deo(
        list(rep(list(md$tree_init$edge), length(temperatures))),
        md$logP, md$logA, niter, 91L, temperatures, 5L)

    # The endpoint is a swap boundary, so this checks the post-swap overwrite
    # rather than merely comparing two states after an update-only iteration.
    expect_gt(deo$swap_accepts[1L, 1L], 0L)
    expect_identical(
        as.integer(deo$traces[[1L]][[niter + 1L]]),
        as.integer(deo$final_states[[1L]][[1L]]))
})

test_that("DEO RNG streams are deterministic across worker counts", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.1, 1.2, 1.3, 1.4)
    starts <- rep(
        list(rep(list(md$tree_init$edge), length(temperatures))),
        3L)
    seeds <- c(11L, 12L, 13L)
    on.exit(RcppParallel::setThreadOptions(numThreads = 2L), add = TRUE)

    RcppParallel::setThreadOptions(numThreads = 1L)
    one_worker <- mitodrift:::tree_mc3_parallel_seeded_deo(
        starts, md$logP, md$logA, rep(100L, 3L), seeds,
        temperatures, 5L)
    RcppParallel::setThreadOptions(numThreads = 15L)
    many_workers <- mitodrift:::tree_mc3_parallel_seeded_deo(
        starts, md$logP, md$logA, rep(100L, 3L), seeds,
        temperatures, 5L)

    expect_identical(many_workers, one_worker)
})

test_that("batched MC3 defaults exactly to the legacy RNN schedule", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.1, 1.2)

    run_case <- function(explicit_scheme) {
        outdir <- tempfile("mitodrift-mc3-rnn-default-")
        dir.create(outdir)
        args <- list(
            phy_init = md$tree_init,
            logP_list = md$logP,
            logA_vec = md$logA,
            outfile = file.path(outdir, "trace.qs2"),
            max_iter = 20L,
            nchains = 2L,
            ncores = 2L,
            batch_size = 10L,
            mc3_temperatures = temperatures,
            mc3_swap_interval = 5L,
            mc3_statefile = file.path(outdir, "state.rds")
        )
        if (explicit_scheme) args$mc3_swap_scheme <- "rnn"
        trace <- do.call(mitodrift:::run_tree_mcmc_batch, args)
        list(trace = trace, checkpoint = readRDS(args$mc3_statefile))
    }

    default <- run_case(FALSE)
    explicit <- run_case(TRUE)

    expect_identical(default, explicit)
    expect_identical(default$checkpoint$swap_scheme, "rnn")
})

test_that("batched DEO resume preserves phase and is bit-identical", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.05, 1.1, 1.15)

    full_dir <- tempfile("mitodrift-mc3-deo-full-")
    resumed_dir <- tempfile("mitodrift-mc3-deo-resumed-")
    dir.create(full_dir)
    dir.create(resumed_dir)

    run_deo <- function(outdir, max_iter, resume = FALSE) {
        trace_file <- file.path(outdir, "trace.qs2")
        state_file <- file.path(outdir, "state.rds")
        trace <- mitodrift:::run_tree_mcmc_batch(
            md$tree_init, md$logP, md$logA, trace_file,
            max_iter = max_iter, nchains = 2L, ncores = 2L,
            batch_size = 5L, resume = resume,
            mc3_temperatures = temperatures,
            mc3_swap_interval = 5L, mc3_swap_scheme = "deo",
            mc3_statefile = state_file)
        list(trace = trace, checkpoint = readRDS(state_file))
    }

    uninterrupted <- run_deo(full_dir, 20L)
    run_deo(resumed_dir, 5L)
    resumed <- run_deo(resumed_dir, 20L, resume = TRUE)

    expect_identical(resumed, uninterrupted)
    expect_identical(resumed$checkpoint$swap_scheme, "deo")
    expect_identical(resumed$checkpoint$deo_phase, 0L)
    expect_identical(
        resumed$checkpoint$swap_attempts,
        matrix(rep(c(2L, 2L, 2L), 2L), nrow = 2L, byrow = TRUE))
})

test_that("MC3 resume rejects a changed swap scheme", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-mc3-scheme-mismatch-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    state_file <- file.path(outdir, "state.rds")
    temperatures <- c(1, 1.1, 1.2)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        max_iter = 5L, nchains = 2L, ncores = 2L,
        batch_size = 5L, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_swap_scheme = "rnn",
        mc3_statefile = state_file, return_trace = FALSE)

    expect_error(
        mitodrift:::run_tree_mcmc_batch(
            md$tree_init, md$logP, md$logA, trace_file,
            max_iter = 10L, nchains = 2L, ncores = 2L,
            batch_size = 5L, resume = TRUE,
            mc3_temperatures = temperatures,
            mc3_swap_interval = 5L, mc3_swap_scheme = "deo",
            mc3_statefile = state_file, return_trace = FALSE),
        "swap scheme does not match"
    )
})

test_that("legacy MC3 checkpoints without a scheme resume as RNN", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-mc3-legacy-scheme-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    state_file <- file.path(outdir, "state.rds")
    temperatures <- c(1, 1.1, 1.2)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        max_iter = 5L, nchains = 2L, ncores = 2L,
        batch_size = 5L, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)
    checkpoint <- readRDS(state_file)
    checkpoint$swap_scheme <- NULL
    saveRDS(checkpoint, state_file)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        max_iter = 10L, nchains = 2L, ncores = 2L,
        batch_size = 5L, resume = TRUE,
        mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)

    expect_identical(readRDS(state_file)$swap_scheme, "rnn")
})
