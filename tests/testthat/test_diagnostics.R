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
