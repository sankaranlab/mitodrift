#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages({
    library(igraph)
    library(dplyr)
    library(ggplot2)
    library(ape)
    library(tidygraph)
    library(data.table)
    library(stringr)
    library(parallel)
    library(tidytree)
    library(glue)
    library(patchwork)
    library(phangorn)
    library(ggtree)
    library(ggraph)
    library(mitodrift)
    library(optparse)
})

message("=== mitodrift package version ===")
message("mitodrift version: ", as.character(utils::packageVersion("mitodrift")), "")

# repo_dir = "/lab-share/Hem-Sankaran-e2/Public/projects/tgao/tools/mitodrift"
# devtools::load_all(repo_dir)

# Define command line options
option_list <- list(
    make_option(
        c("-m", "--mut_dat"),
        type = "character",
        help = "Mutation data file (CSV format with columns: variant, cell, d, a)",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-A", "--amat"),
        type = "character",
        help = "Alternative allele count matrix (CSV; first column = variant)",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-D", "--dmat"),
        type = "character",
        help = "Total depth matrix (CSV; first column = variant)",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        default = "mitodrift_results",
        help = "Output directory for results",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-p", "--ncores"),
        type = "integer",
        default = 1,
        help = "Number of cores to use",
        metavar = "INTEGER"
    ),
    make_option(
        c("--ncores_nj"),
        type = "integer",
        default = 1,
        help = "Number of cores to use for initial NJ tree construction via ParDist; notice this may change the resulting tree",
        metavar = "INTEGER"
    ),
    make_option(
        c("-u", "--ncores_em"),
        type = "integer",
        default = 1,
        help = "Number of cores to use for EM parameter fitting",
        metavar = "INTEGER"
    ),
    make_option(
        c("-q", "--ncores_qs"),
        type = "integer",
        default = 1,
        help = "Number of cores to use for file saving/writing (qs2) in MCMC and tree annotation",
        metavar = "INTEGER"
    ),
    make_option(
        c("--ncores_annot"),
        type = "integer",
        default = NA_integer_,
        help = "Number of cores to use for clade annotation (default: same as --ncores)",
        metavar = "INTEGER"
    ),
    make_option(
        c("-k", "--k"),
        type = "integer",
        default = 20,
        help = "Number of VAF bins",
        metavar = "INTEGER"
    ),
    make_option(
        c("-n", "--npop"),
        type = "integer",
        default = 600,
        help = "Population size",
        metavar = "INTEGER"
    ),
    make_option(
        c("-e", "--eps"),
        type = "double",
        default = 0.001,
        help = "Mutation rate per branch; automatically fitted using EM if fit_params=TRUE",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-s", "--err"),
        type = "double",
        default = 0,
        help = "Variant detection error rate; automatically fitted using EM if fit_params=TRUE",
        metavar = "DOUBLE"
    ),
    make_option(
        c("--contamination_rate"),
        type = "double",
        default = 0,
        help = "Global mitochondrial transfer/cross-contamination fraction (default: 0)",
        metavar = "DOUBLE"
    ),
    make_option(
        c("--fit_contamination"),
        type = "logical",
        default = FALSE,
        help = "Whether to estimate contamination_rate during EM (default: FALSE)",
        metavar = "LOGICAL"
    ),
    make_option(
        c("--contamination_max"),
        type = "double",
        default = 0.2,
        help = "Upper bound for contamination_rate during EM (default: 0.2)",
        metavar = "DOUBLE"
    ),
    make_option(
        c("--contamination_ecm_steps"),
        type = "integer",
        default = 1,
        help = "Conditional err/contamination updates per EM iteration (default: 1)",
        metavar = "INTEGER"
    ),
    make_option(
        c("-g", "--ngen"),
        type = "integer",
        default = 100,
        help = "Number of generations; automatically fitted using EM if fit_params=TRUE",
        metavar = "INTEGER"
    ),
    make_option(
        c("-f", "--fit_params"),
        type = "logical",
        default = TRUE,
        help = "Whether to fit parameters using EM",
        metavar = "LOGICAL"
    ),
    make_option(
        c("-t", "--fit_param_max_iter"),
        type = "integer",
        default = 100,
        help = "Maximum EM iterations for parameter fitting",
        metavar = "INTEGER"
    ),
    make_option(
        c("-c", "--fit_param_epsilon"),
        type = "double",
        default = 1e-3,
        help = "Convergence threshold for EM parameter fitting",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-i", "--ml_iter"),
        type = "integer",
        default = 0,
        help = "Maximum iterations for tree optimization",
        metavar = "INTEGER"
    ),
    make_option(
        c("-j", "--tree_mcmc_iter"),
        type = "integer",
        default = 10000,
        help = "Maximum iterations for phylogenetic MCMC",
        metavar = "INTEGER"
    ),
    make_option(
        c("-l", "--tree_mcmc_chains"),
        type = "integer",
        default = 1,
        help = "Number of MCMC chains for phylogenetic sampling",
        metavar = "INTEGER"
    ),
    make_option(
        c("-b", "--tree_mcmc_burnin"),
        type = "integer",
        default = 0,
        help = "Burnin for phylogenetic MCMC",
        metavar = "INTEGER"
    ),
    make_option(
        c("-d", "--tree_mcmc_batch_size"),
        type = "integer",
        default = 1000,
        help = "Batch size for phylogenetic MCMC",
        metavar = "INTEGER"
    ),
    make_option(
        c("-y", "--tree_mcmc_diag"),
        type = "logical",
        default = TRUE,
        help = "Whether to run diagnostics (e.g., ASDSF) after each MCMC batch",
        metavar = "LOGICAL"
    ),
    make_option(
        c("--tree_mcmc_checkpoint_every"),
        type = "integer",
        default = 1L,
        help = "Persist MCMC trace/state every N diagnostic batches (final batch always saved)",
        metavar = "INTEGER"
    ),
    make_option(
        c("-z", "--conv_thres"),
        type = "double",
        default = NULL,
        help = "Convergence threshold for ASDSF",
        metavar = "DOUBLE"
    ),
    make_option(
        c("--mc3_temperatures"),
        type = "character",
        default = NULL,
        help = "Comma-separated MC3 temperature ladder starting at 1; omitted uses independent MCMC",
        metavar = "CHARACTER"
    ),
    make_option(
        c("--mc3_swap_interval"),
        type = "integer",
        default = 10L,
        help = "Iterations between MC3 adjacent-temperature swap proposals",
        metavar = "INTEGER"
    ),
    make_option(
        c("--mc3_ncores"),
        type = "integer",
        default = NULL,
        help = "MC3 sampling threads; with allocated CPUs, use up to chains times temperatures",
        metavar = "INTEGER"
    ),
    make_option(
        c("-r", "--resume"),
        type = "logical",
        default = FALSE,
        help = "Whether to resume from existing files",
        metavar = "LOGICAL"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

if (!is.null(opts$mc3_temperatures)) {
    opts$mc3_temperatures <- as.numeric(strsplit(opts$mc3_temperatures, ',', fixed = TRUE)[[1]])
    if (anyNA(opts$mc3_temperatures)) stop('--mc3_temperatures must be comma-separated numbers')
}

# Default annotation cores to ncores if not set.
if (is.na(opts$ncores_annot)) {
    opts$ncores_annot <- opts$ncores
}

# Print parameters
message("=== MitoDrift Analysis Parameters ===")
for (arg in names(opts)) {
    value <- opts[[arg]]
    message(arg, ": ", if (length(value) > 1L) paste(value, collapse = ",") else value)
}
if (opts$resume) {
    message("RESUME MODE: Will attempt to resume from existing files")
}

# Check required parameters
if (!is.null(opts$mut_dat)) {
    if (!is.null(opts$amat) || !is.null(opts$dmat)) {
        stop("Provide either --mut_dat or both --amat and --dmat, not both")
    }
} else {
    if (is.null(opts$amat) || is.null(opts$dmat)) {
        stop("Must provide --mut_dat or both --amat and --dmat")
    }
}

# Create output directory
if (!dir.exists(opts$outdir)) {
    dir.create(opts$outdir, recursive = TRUE)
    message("Created output directory: ", opts$outdir)
}

# Define output file paths
ml_trace_file <- file.path(opts$outdir, "tree_ml_trace.qs2")
mcmc_trace_file <- file.path(opts$outdir, "tree_mcmc_trace.qs2")
mcmc_diag_file <- file.path(opts$outdir, "tree_mcmc_diag.rds")
annot_tree_file <- file.path(opts$outdir, "tree_annotated.newick")
mitodrift_object_file <- file.path(opts$outdir, "mitodrift_object.rds")

start_time <- Sys.time()

# Step 1: Load mutation data / matrices
message("\n=== Creating MitoDrift object ===")

load_matrix_file <- function(path) {
    dt <- fread(path)
    if (ncol(dt) < 2) {
        stop(glue("Matrix file '{path}' must include a row identifier column plus data columns"))
    }
    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- dt[[1]]
    storage.mode(mat) <- "numeric"
    mat
}

if (!is.null(opts$mut_dat)) {
    mut_dat <- fread(opts$mut_dat) %>% 
        group_by(variant) %>%
        filter(sum(a > 0) > 1) %>%
        ungroup()
    message("Input mode: long-format mutation table")
    message("Number of unique variants: ", length(unique(mut_dat$variant)))
    message("Number of unique cells: ", length(unique(mut_dat$cell)))
} else {
    amat <- load_matrix_file(opts$amat)
    dmat <- load_matrix_file(opts$dmat)
    
    if (!identical(dim(amat), dim(dmat))) {
        stop("amat and dmat must have identical dimensions")
    }
    # Check for NA values in input matrices
    if (any(is.na(amat)) || any(is.na(dmat))) {
        stop('Input matrices cannot contain NA values')
    }

    keep_variants <- rowSums(amat > 0) > 1
    if (!any(keep_variants)) {
        stop("No variants remain in amat after filtering for >1 non-zero counts")
    }
    n_removed <- nrow(amat) - sum(keep_variants)
    if (n_removed > 0) {
        message("Filtered ", n_removed, " variants with <= 1 non-zero count")
    }
    amat <- amat[keep_variants, , drop = FALSE]
    dmat <- dmat[keep_variants, , drop = FALSE]
    message("Input mode: matrix pair")
    message("Matrix dimensions: ", nrow(amat), " variants x ", ncol(amat), " cells")
}

# Step 2: Create MitoDrift object
if (opts$resume) {
    message("=== Resuming from existing MitoDrift object ===")
    saved_md <- readRDS(mitodrift_object_file)
    if (is.null(saved_md$amat) || is.null(saved_md$dmat)) {
        stop("Saved MitoDrift object lacks `amat`/`dmat`; cannot resume.")
    }
    md <- MitoDrift$new(
        amat = saved_md$amat,
        dmat = saved_md$dmat,
        model_params = saved_md$model_params,
        build_tree = FALSE,
        ncores = opts$ncores_nj
    )
    md <- md$copy(saved_md)
} else {
    model_params <- c(
        eps = opts$eps,
        err = opts$err,
        contamination_rate = opts$contamination_rate,
        npop = opts$npop,
        ngen = opts$ngen,
        k = opts$k
    )

    if (!is.null(opts$mut_dat)) {
        message(glue("Building NJ tree for long-format input using {opts$ncores_nj} core(s)"))
        md <- MitoDrift$new(
            mut_dat = mut_dat,
            model_params = model_params,
            ncores = opts$ncores_nj
        )
    } else {
        message(glue("Building NJ tree for matrix input using {opts$ncores_nj} core(s)"))
        md <- MitoDrift$new(
            amat = amat,
            dmat = dmat,
            model_params = model_params,
            ncores = opts$ncores_nj
        )
    }
    saveRDS(md, mitodrift_object_file)
}

if (!(!is.null(md$logP) && !is.null(md$logA))) {
    if (opts$fit_params) {
        message("\n=== Fitting parameters using EM ===")
        fitted_params <- md$fit_params_em(
            tree_fit = md$tree_init,
            max_iter = opts$fit_param_max_iter,
            epsilon = opts$fit_param_epsilon,
            ncores = opts$ncores_em,
            fit_contamination = opts$fit_contamination,
            contamination_bounds = c(0, opts$contamination_max),
            emission_ecm_steps = opts$contamination_ecm_steps
        )
    } else {
        message("\n=== Skipping parameter fitting ===")
    }

    message("\n=== Initializing model components ===")
    md$make_model(ncores = opts$ncores)
    saveRDS(md, mitodrift_object_file)
} else {
    message("\n=== Model components already initialized; skipping parameter fitting and initialization ===")
}

message("\n=== Optimizing tree ===")
md$optimize_tree(
    max_iter = opts$ml_iter,
    ncores = opts$ncores,
    outfile = ml_trace_file,
    resume = opts$resume
)

saveRDS(md, mitodrift_object_file)

message("\n=== Running phylogenetic MCMC ===")
md$run_mcmc(
    max_iter = opts$tree_mcmc_iter,
    nchains = opts$tree_mcmc_chains,
    ncores = opts$ncores,
    ncores_qs = opts$ncores_qs,
    batch_size = opts$tree_mcmc_batch_size,
    diag = opts$tree_mcmc_diag,
    conv_thres = opts$conv_thres,
    outfile = mcmc_trace_file,
    diagfile = mcmc_diag_file,
    resume = opts$resume,
    mc3_temperatures = opts$mc3_temperatures,
    mc3_swap_interval = opts$mc3_swap_interval,
    mc3_ncores = opts$mc3_ncores,
    mc3_statefile = file.path(opts$outdir, "tree_mc3_state.rds"),
    checkpoint_every = opts$tree_mcmc_checkpoint_every
)

saveRDS(md, mitodrift_object_file)

message("\n=== Annotating tree with clade frequencies ===")
md$annotate_tree(
    burnin = opts$tree_mcmc_burnin,
    ncores = opts$ncores_annot,
    ncores_qs = opts$ncores_qs
)

write.tree(md$tree_annot, annot_tree_file)
saveRDS(md, mitodrift_object_file)

message("\n=== MitoDrift analysis completed successfully! ===")

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "secs")
message(glue("Total runtime: {format(elapsed)} (start: {start_time}, end: {end_time})"))
