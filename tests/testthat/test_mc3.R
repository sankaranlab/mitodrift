test_that("one-temperature MC3 exactly matches the original sampler", {
    md <- make_mc3_test_model()
    niter <- 100L
    seed <- 42L

    original <- mitodrift:::tree_mcmc_cpp_cached_threadsafe(
        md$tree_init$edge, md$logP, md$logA, niter, seed, TRUE)
    coupled <- mitodrift:::tree_mc3_parallel_seeded(
        list(list(md$tree_init$edge)), md$logP, md$logA,
        niter, seed, 1, 10L)

    expect_identical(coupled$traces[[1]], original)
    expect_equal(dim(coupled$swap_attempts), c(1L, 0L))
})

test_that("MC3 records adjacent swaps and returns valid cold and heated states", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.5, 2.25, 3.375)
    niter <- 200L
    swap_interval <- 5L

    coupled <- mitodrift:::tree_mc3_parallel_seeded(
        list(rep(list(md$tree_init$edge), length(temperatures))),
        md$logP, md$logA, niter, 7L, temperatures, swap_interval)

    expect_length(coupled$traces[[1]], niter + 1L)
    expect_length(coupled$final_states[[1]], length(temperatures))
    expect_equal(sum(coupled$swap_attempts), niter %/% swap_interval)
    expect_true(all(coupled$swap_accepts >= 0L))
    expect_true(all(coupled$swap_accepts <= coupled$swap_attempts))

    scores <- vapply(coupled$final_states[[1]], function(edges) {
        mitodrift:::score_tree_bp_wrapper2(edges, md$logP, md$logA)
    }, numeric(1))
    expect_true(all(is.finite(scores)))
})

test_that("MC3 cold samples match the exact posterior over five-tip topologies", {
    md <- make_mc3_test_model(n_tips = 5L)
    exact <- exact_rooted_topology_posterior(md)
    expect_length(exact$trees, 105L)
    expect_equal(sum(exact$probability), 1, tolerance = 1e-14)

    temperatures <- c(1, 1.5, 2.25)
    nchains <- 4L
    niter <- 20000L
    burnin <- 4000L
    starts <- rep(
        list(rep(list(md$tree_init$edge), length(temperatures))),
        nchains)
    coupled <- mitodrift:::tree_mc3_parallel_seeded(
        starts, md$logP, md$logA,
        rep(niter, nchains), seq_len(nchains), temperatures, 10L)
    estimate <- estimate_topology_posterior(
        coupled$traces, exact$keys, 5L, burnin)

    total_variation <- 0.5 * sum(abs(estimate - exact$probability))
    max_abs_error <- max(abs(estimate - exact$probability))
    expect_lt(total_variation, 0.02)
    expect_lt(max_abs_error, 0.01)
    expect_gt(stats::cor(estimate, exact$probability), 0.999)
})

test_that("flattened MC3 RNG streams are deterministic across thread counts", {
    md <- make_mc3_test_model()
    temperatures <- c(1, 1.5, 2.25, 3.375)
    starts <- rep(list(rep(list(md$tree_init$edge), length(temperatures))), 3L)
    seeds <- c(11L, 12L, 13L)
    on.exit(RcppParallel::setThreadOptions(numThreads = 2L), add = TRUE)

    RcppParallel::setThreadOptions(numThreads = 1L)
    serial_schedule <- mitodrift:::tree_mc3_parallel_seeded(
        starts, md$logP, md$logA, rep(100L, 3L), seeds,
        temperatures, 5L)
    RcppParallel::setThreadOptions(numThreads = 12L)
    parallel_schedule <- mitodrift:::tree_mc3_parallel_seeded(
        starts, md$logP, md$logA, rep(100L, 3L), seeds,
        temperatures, 5L)

    expect_identical(parallel_schedule, serial_schedule)
})

test_that("MC3 validates its temperature ladder", {
    md <- make_mc3_test_model()
    starts <- list(list(md$tree_init$edge, md$tree_init$edge))

    expect_error(
        mitodrift:::tree_mc3_parallel_seeded(starts, md$logP, md$logA, 10L, 1L, c(1.1, 2), 5L),
        "start at 1"
    )
    expect_error(
        mitodrift:::tree_mc3_parallel_seeded(starts, md$logP, md$logA, 10L, 1L, c(1, 1), 5L),
        "strictly increase"
    )
})

test_that("batched MC3 checkpoints and resumes heated states", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-mc3-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    state_file <- file.path(outdir, "state.rds")
    diag_file <- file.path(outdir, "diag.rds")

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 20L, nchains = 2L, ncores = 2L,
        batch_size = 10L, mc3_temperatures = c(1, 1.5, 2.25),
        mc3_swap_interval = 5L, mc3_statefile = state_file)
    first_checkpoint <- readRDS(state_file)
    expect_equal(first_checkpoint$completed_iters, 20L)
    expect_equal(sum(first_checkpoint$swap_attempts), 8L)

    mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA, trace_file,
        diagfile = diag_file, max_iter = 30L, nchains = 2L, ncores = 2L,
        batch_size = 10L, resume = TRUE, mc3_temperatures = c(1, 1.5, 2.25),
        mc3_swap_interval = 5L, mc3_statefile = state_file)
    resumed_trace <- qs2::qd_read(trace_file)
    resumed_checkpoint <- readRDS(state_file)
    expect_equal(lengths(resumed_trace), c(`1` = 31L, `2` = 31L))
    expect_equal(resumed_checkpoint$completed_iters, 30L)
    expect_equal(sum(resumed_checkpoint$swap_attempts), 12L)
})

test_that("a one-temperature batch run is exactly ordinary MCMC", {
    md <- make_mc3_test_model()
    ordinary_dir <- tempfile("mitodrift-ordinary-")
    one_temp_dir <- tempfile("mitodrift-one-temp-")
    dir.create(ordinary_dir)
    dir.create(one_temp_dir)

    ordinary <- mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA,
        outfile = file.path(ordinary_dir, "trace.qs2"),
        max_iter = 20L, nchains = 2L, ncores = 2L, batch_size = 10L)
    one_temp <- mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA,
        outfile = file.path(one_temp_dir, "trace.qs2"),
        max_iter = 20L, nchains = 2L, ncores = 2L, batch_size = 10L,
        mc3_temperatures = 1,
        mc3_statefile = file.path(one_temp_dir, "state.rds"))

    expect_identical(one_temp, ordinary)
    expect_false(file.exists(file.path(one_temp_dir, "state.rds")))
})

test_that("deferred checkpoints retain the complete final trace and diagnostics", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-checkpoint-cadence-")
    dir.create(outdir)
    trace_file <- file.path(outdir, "trace.qs2")
    diag_file <- file.path(outdir, "diag.rds")

    result <- mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA,
        outfile = trace_file, diagfile = diag_file,
        max_iter = 25L, nchains = 2L, ncores = 2L,
        batch_size = 10L, checkpoint_every = 2L)

    persisted <- qs2::qd_read(trace_file)
    diagnostics <- readRDS(diag_file)
    expect_equal(lengths(result), c(`1` = 26L, `2` = 26L))
    expect_identical(persisted, result)
    expect_equal(diagnostics$completed_iters, c(10L, 20L, 25L))
})
