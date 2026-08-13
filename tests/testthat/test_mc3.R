make_mc3_test_model <- function() {
    mut_dat <- utils::read.csv(system.file("extdata", "small_test_mut_dat.csv", package = "mitodrift"))
    md <- MitoDrift$new(
        mut_dat = mut_dat,
        model_params = c(eps = 0.001, err = 0, npop = 600, ngen = 100, k = 20),
        ncores = 1
    )
    md$make_model(ncores = 1)
    md
}

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
