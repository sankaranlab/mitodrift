test_that("incremental ASDSF matches a full-history recomputation", {
    md <- make_mc3_test_model()
    trees <- mitodrift:::tree_mcmc_parallel_seeded(
        rep(list(md$tree_init$edge), 2), md$logP, md$logA,
        c(20L, 20L), c(1L, 2L))
    trees <- lapply(trees, mitodrift:::restore_elist)

    first <- lapply(trees, function(chain) chain[1:11])
    second <- lapply(trees, function(chain) chain[12:21])
    state <- mitodrift:::initialize_target_tree_asdsf_state(md$tree_init, first)
    state <- mitodrift:::update_target_tree_asdsf_state(state, md$tree_init, second)

    expect_equal(
        mitodrift:::target_tree_asdsf_from_state(state),
        mitodrift:::compute_target_tree_asdsf(md$tree_init, trees),
        tolerance = 0
    )
})

test_that("convergence mode is not capped by max_iter", {
    md <- make_mc3_test_model()
    outdir <- tempfile("mitodrift-convergence-cap-")
    dir.create(outdir)

    result <- mitodrift:::run_tree_mcmc_batch(
        md$tree_init, md$logP, md$logA,
        outfile = file.path(outdir, "trace.qs2"),
        diagfile = file.path(outdir, "diag.rds"),
        max_iter = 1L, conv_thres = 1, nchains = 2L, ncores = 2L,
        batch_size = 10L)

    expect_equal(lengths(result), c(`1` = 11L, `2` = 11L))
    expect_equal(readRDS(file.path(outdir, "diag.rds"))$completed_iters, 10L)
})
