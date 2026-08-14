#' Compute depth-weighted global VAFs
#'
#' Computes one pooled alternate-allele fraction per variant for use as the
#' source-pool VAF in the contamination-aware emission model. Variants with no
#' observed depth are assigned a global VAF of zero because they contribute no
#' read likelihood.
#'
#' @param amat Alternative-allele count matrix (variants by cells).
#' @param dmat Total-depth matrix with the same dimensions as `amat`.
#' @return A named numeric vector with one global VAF per variant.
#' @keywords internal
.compute_global_vaf <- function(amat, dmat) {
    if (!is.matrix(amat) || !is.matrix(dmat) || !identical(dim(amat), dim(dmat))) {
        stop("amat and dmat must be matrices with identical dimensions")
    }
    if (any(!is.finite(amat)) || any(!is.finite(dmat)) ||
        any(amat < 0) || any(dmat < 0) || any(amat > dmat)) {
        stop("Read counts must be finite and satisfy 0 <= amat <= dmat")
    }

    alt_sum <- rowSums(amat)
    depth_sum <- rowSums(dmat)
    global_vaf <- numeric(length(depth_sum))
    covered <- depth_sum > 0
    global_vaf[covered] <- alt_sum[covered] / depth_sum[covered]
    names(global_vaf) <- rownames(amat)
    global_vaf
}

#' Validate a scalar contamination rate
#'
#' @param contamination_rate Fraction of molecules or reads drawn from the
#'   global contamination pool.
#' @return The validated scalar rate.
#' @keywords internal
.validate_contamination_rate <- function(contamination_rate) {
    if (length(contamination_rate) != 1L || !is.finite(contamination_rate) ||
        contamination_rate < 0 || contamination_rate > 1) {
        stop("contamination_rate must be one finite value in [0, 1]")
    }
    as.numeric(contamination_rate)
}
