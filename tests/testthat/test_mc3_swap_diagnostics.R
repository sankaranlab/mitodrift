test_that("MC3 diagnostics retain exact batch and adjacent-pair swap counts", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-mc3-swap-diagnostics-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    state_file <- file.path(outdir, "state.rds")
    diag_file <- file.path(outdir, "diag.rds")
    temperatures <- c(1, 1.25, 1.75)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 20L, nchains = 2L, ncores = 2L,
        batch_size = 10L, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)

    diagnostics <- readRDS(diag_file)
    checkpoint <- readRDS(state_file)

    expect_equal(nrow(diagnostics), 2L)
    expect_equal(diagnostics$mc3_swap_attempts_batch, c(4L, 4L))
    expect_equal(
        diagnostics$mc3_swap_acceptance_batch,
        diagnostics$mc3_swap_accepts_batch / diagnostics$mc3_swap_attempts_batch)
    for (batch_id in seq_len(nrow(diagnostics))) {
        pair_attempts <- diagnostics$mc3_swap_attempts_by_pair_batch[[batch_id]]
        pair_accepts <- diagnostics$mc3_swap_accepts_by_pair_batch[[batch_id]]
        expect_equal(
            diagnostics$mc3_swap_acceptance_by_pair_batch[[batch_id]],
            ifelse(pair_attempts > 0L, pair_accepts / pair_attempts, NA_real_))
    }
    expect_equal(
        unname(Reduce(`+`, diagnostics$mc3_swap_attempts_by_pair_batch)),
        colSums(checkpoint$swap_attempts))
    expect_equal(
        unname(Reduce(`+`, diagnostics$mc3_swap_accepts_by_pair_batch)),
        colSums(checkpoint$swap_accepts))
    expect_true(all(vapply(
        diagnostics$mc3_swap_attempts_by_pair_batch,
        function(x) identical(names(x), c("pair_1", "pair_2")),
        logical(1))))
    expected_pairs <- matrix(
        c(1, 1.25, 1.25, 1.75), ncol = 2L,
        dimnames = list(c("pair_1", "pair_2"), c("lower", "upper")))
    expect_true(all(vapply(
        diagnostics$mc3_swap_pair_temperatures,
        function(x) isTRUE(all.equal(x, expected_pairs)),
        logical(1))))

    cumulative_acceptance <- sum(checkpoint$swap_accepts) / sum(checkpoint$swap_attempts)
    expect_equal(tail(diagnostics$mc3_swap_acceptance, 1L), cumulative_acceptance)
})

test_that("resumed MC3 diagnostics record only newly attempted swaps", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-mc3-swap-diagnostics-resume-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    state_file <- file.path(outdir, "state.rds")
    diag_file <- file.path(outdir, "diag.rds")
    temperatures <- c(1, 1.25, 1.75)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 10L, nchains = 2L, ncores = 2L,
        batch_size = 10L, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)
    before <- readRDS(state_file)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 20L, nchains = 2L, ncores = 2L,
        batch_size = 10L, resume = TRUE, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)
    diagnostics <- readRDS(diag_file)
    after <- readRDS(state_file)

    expect_equal(diagnostics$completed_iters, c(10L, 20L))
    expect_equal(diagnostics$mc3_swap_attempts_batch, c(4L, 4L))
    expect_equal(
        unname(diagnostics$mc3_swap_attempts_by_pair_batch[[2L]]),
        colSums(after$swap_attempts - before$swap_attempts))
    expect_equal(
        unname(diagnostics$mc3_swap_accepts_by_pair_batch[[2L]]),
        colSums(after$swap_accepts - before$swap_accepts))
    expect_equal(
        diagnostics$mc3_swap_acceptance[[1L]],
        sum(before$swap_accepts) / sum(before$swap_attempts))
    expect_equal(
        diagnostics$mc3_swap_acceptance[[2L]],
        sum(after$swap_accepts) / sum(after$swap_attempts))
})

test_that("MC3 resume accepts diagnostics written before batch swap columns", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-mc3-legacy-swap-diagnostics-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    state_file <- file.path(outdir, "state.rds")
    diag_file <- file.path(outdir, "diag.rds")
    temperatures <- c(1, 1.25, 1.75)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 10L, nchains = 2L, ncores = 2L,
        batch_size = 10L, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)
    legacy_diagnostics <- readRDS(diag_file)[
        , c("batch", "completed_iters", "asdsf", "mc3_swap_acceptance")]
    saveRDS(legacy_diagnostics, diag_file)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 20L, nchains = 2L, ncores = 2L,
        batch_size = 10L, resume = TRUE, mc3_temperatures = temperatures,
        mc3_swap_interval = 5L, mc3_statefile = state_file,
        return_trace = FALSE)
    diagnostics <- readRDS(diag_file)

    expect_equal(diagnostics$completed_iters, c(10L, 20L))
    expect_true(is.na(diagnostics$mc3_swap_attempts_batch[[1L]]))
    expect_equal(diagnostics$mc3_swap_attempts_batch[[2L]], 4L)
    expect_null(diagnostics$mc3_swap_attempts_by_pair_batch[[1L]])
    expect_equal(sum(diagnostics$mc3_swap_attempts_by_pair_batch[[2L]]), 4L)
})
