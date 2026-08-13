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
