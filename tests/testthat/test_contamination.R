make_contamination_test_counts <- function() {
    amat <- matrix(
        as.integer(c(0, 2, 5, 1, 4, 8)),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("v1", "v2"), c("c1", "c2", "c3"))
    )
    dmat <- matrix(
        as.integer(c(10, 10, 10, 10, 10, 10)),
        nrow = 2,
        byrow = TRUE,
        dimnames = dimnames(amat)
    )
    list(amat = amat, dmat = dmat)
}

test_that("global VAF is pooled by read depth", {
    counts <- make_contamination_test_counts()
    expect_equal(
        .compute_global_vaf(counts$amat, counts$dmat),
        c(v1 = 7 / 30, v2 = 13 / 30)
    )

    invalid <- counts$amat
    invalid[1, 1] <- 11L
    expect_error(
        .compute_global_vaf(invalid, counts$dmat),
        "0 <= amat <= dmat",
        fixed = TRUE
    )
})

test_that("zero contamination exactly preserves leaf likelihoods", {
    counts <- make_contamination_test_counts()
    vafs <- c(0, 0.5, 1)
    global_vaf <- .compute_global_vaf(counts$amat, counts$dmat)

    original <- get_leaf_liks_mat_cpp(
        counts$amat, counts$dmat, vafs, eps = 0.01, log = TRUE
    )
    extended <- get_leaf_liks_mat_cpp(
        counts$amat, counts$dmat, vafs, eps = 0.01, log = TRUE,
        contamination_rate = 0, global_vaf = global_vaf
    )

    expect_identical(extended, original)
})

test_that("contamination mixes latent and global VAF before read error", {
    counts <- make_contamination_test_counts()
    vafs <- c(0, 0.5, 1)
    err <- 0.01
    contamination_rate <- 0.2
    global_vaf <- .compute_global_vaf(counts$amat, counts$dmat)

    observed <- get_leaf_liks_mat_cpp(
        counts$amat, counts$dmat, vafs, eps = err, log = TRUE,
        contamination_rate = contamination_rate, global_vaf = global_vaf
    )

    for (variant in seq_len(nrow(counts$amat))) {
        probability <- pmin(
            (1 - contamination_rate) * vafs +
                contamination_rate * global_vaf[[variant]] + err,
            1 - err
        )
        expected <- vapply(seq_len(ncol(counts$amat)), function(cell) {
            dbinom(
                counts$amat[variant, cell],
                counts$dmat[variant, cell],
                probability,
                log = TRUE
            )
        }, numeric(length(vafs)))
        expect_equal(unname(observed[[variant]]), unname(expected), tolerance = 1e-14)
    }

    expect_error(
        get_leaf_liks_mat_cpp(
            counts$amat, counts$dmat, vafs,
            contamination_rate = contamination_rate
        ),
        "global_vaf is required"
    )
})

test_that("MitoDrift defaults contamination to zero", {
    counts <- make_contamination_test_counts()
    model_params <- c(eps = 0.001, err = 0.01, npop = 20, ngen = 2, k = 2)
    model <- MitoDrift$new(
        amat = counts$amat,
        dmat = counts$dmat,
        model_params = model_params,
        build_tree = FALSE
    )
    model$tree_init <- ape::read.tree(text = "((c1,c2),c3);")
    model$make_model()

    expect_equal(model$model_params[["contamination_rate"]], 0)
    expect_equal(model$global_vaf, .compute_global_vaf(counts$amat, counts$dmat))
})

test_that("EM can conditionally fit one global contamination rate", {
    counts <- make_contamination_test_counts()
    tree <- ape::read.tree(text = "((c1,c2),c3);")

    fit <- suppressMessages(fit_params_em_cpp(
        tree_fit = tree,
        amat = counts$amat,
        dmat = counts$dmat,
        max_iter = 1,
        k = 2,
        npop = 20,
        trace = FALSE,
        fit_contamination = TRUE,
        contamination_bounds = c(0, 0.3)
    ))

    expect_named(fit, c("ngen", "eps", "err", "contamination_rate"))
    expect_true(all(is.finite(fit)))
    expect_gte(fit[["contamination_rate"]], 0)
    expect_lte(fit[["contamination_rate"]], 0.3)
})
