test_that("long and matrix mutation formats round-trip", {
    mut_dat <- data.frame(
        variant = rep(c("v1", "v2"), each = 3),
        cell = rep(c("c1", "c2", "c3"), 2),
        a = c(0, 1, 2, 3, 0, 1),
        d = rep(10, 6)
    )

    amat <- long_to_mat(mut_dat, "a")
    dmat <- long_to_mat(mut_dat, "d")
    round_trip <- mat_to_long(amat, dmat)
    round_trip <- as.data.frame(round_trip)
    round_trip$cell <- as.character(round_trip$cell)

    expected <- mut_dat[order(mut_dat$variant, mut_dat$cell), ]
    observed <- round_trip[order(round_trip$variant, round_trip$cell), ]
    rownames(expected) <- NULL
    rownames(observed) <- NULL
    expect_equal(observed, expected)
})
