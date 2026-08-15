#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(glue)
    library(ape)
    library(dplyr)
    library(ggplot2)
    library(parallel)
    library(mitodrift)
    library(igraph)
    library(CRF)
    if (requireNamespace("qs2", quietly = TRUE)) {
        library(qs2)
    } else {
        stop("Package 'qs2' is required but not installed.")
    }
    if (requireNamespace("digest", quietly = TRUE)) {
        library(digest)
    } else {
        stop("Package 'digest' is required but not installed.")
    }
    if (requireNamespace("matrixStats", quietly = TRUE)) {
        library(matrixStats)
    } else {
        stop("Package 'matrixStats' is required but not installed.")
    }
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: test_mitodrift_mcmc.R <outdir> [ncores]")
}

outdir <- args[1]
ncores <- if (length(args) >= 2) as.integer(args[2]) else 1

# Construct file paths
obj_file <- file.path(outdir, "mitodrift_object.rds")
trace_file <- file.path(outdir, "tree_mcmc_trace.qs2")

if (!file.exists(obj_file)) stop(glue("mitodrift_object.rds not found in {outdir}"))

# Load data
message("Loading mitodrift object...")
md <- readRDS(obj_file)
message("Loading MCMC trace...")
res_batch <- mitodrift:::safe_read_chain(trace_file)
if (is.null(res_batch)) stop(glue("tree_mcmc_trace.qs2 (or its .blocks/ directory) not found in {outdir}"))

A <- md$A
liks <- md$leaf_likelihoods %>% lapply(exp)

# Collect chains
message("Collecting chains...")
mcmc_trees <- mitodrift::collect_chains(res_batch, md$tree_init, burnin = 1000)

# Generate hashes
message("Hashing trees...")
hashes_mcmc <- mclapply(
    mc.cores = ncores,
    mcmc_trees,
    function(tree) {
        digest::digest(write.tree(tree, file = ""), algo = "md5")
    }
) %>% unlist()

hash_table <- table(hashes_mcmc)
unique_trees <- mcmc_trees[match(names(hash_table), hashes_mcmc)]
count_dict <- setNames(as.integer(hash_table), seq_along(unique_trees))

# Compute scores
message("Computing scores for ", length(unique_trees), " unique trees...")
scores <- mclapply(
    mc.cores = ncores,
    seq_along(unique_trees),
    function(i) {
        if (i %% 100 == 0) message("Processing tree ", i)
        
        utree <- unique_trees[[i]]
        
        # Calculate log likelihood
        logL <- mitodrift:::decode_tree(utree, A, liks, score_only = TRUE)$logZ %>% sum()
        
        data.frame(i = i, logL = logL)
    }
) %>% bind_rows()

# Aggregate counts and probabilities
tree_counts <- scores %>% 
    mutate(freq = count_dict[i]) %>%
    arrange(-logL) %>%
    mutate(p = exp(logL - matrixStats::logSumExp(logL))) %>%
    mutate(frac = freq / sum(freq))

# Plotting
message("Generating plot...")
p <- tree_counts %>%
    ggplot(aes(x = p, y = frac)) +
    geom_abline(linetype = 'dashed') +
    geom_point(size = 2) +
    # scale_x_log10() + 
    # scale_y_log10() +
    labs(
        x = "Posterior Probability (calculated)",
        y = "MCMC Frequency",
        title = "MCMC Alignment Plot"
    ) +
    theme_bw()

# Save plot
plot_file <- file.path(outdir, "mcmc_alignment_plot.pdf")
ggsave(plot_file, p, width = 6, height = 5)
message("Plot saved to ", plot_file)

# Save dataframe
df_file <- file.path(outdir, "mcmc_alignment_data.csv")
write.csv(tree_counts, df_file, row.names = FALSE)
message("Dataframe saved to ", df_file)

# Calculate and print metrics
correlation <- cor(tree_counts$p, tree_counts$frac)
r_squared <- correlation^2
rmse <- sqrt(mean((tree_counts$p - tree_counts$frac)^2))

message("\nAgreement Metrics:")
message(sprintf("Correlation (r): %.4f", correlation))
message(sprintf("R-squared: %.4f", r_squared))
message(sprintf("RMSE: %.4f", rmse))
