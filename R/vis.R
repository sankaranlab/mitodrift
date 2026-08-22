#' Plot a phylogenetic tree with VAF heatmap and annotations
#'
#' @param tree A `phylo` object or `tbl_graph` tree.
#' @param df_var Optional data frame of variant calls with columns `cell`, `variant`, and either `vaf`
#'   or `a`/`d` for computing VAF.
#' @param branch_width Numeric branch line width for the tree.
#' @param root_edge Logical; included for compatibility with earlier calls.
#' @param dot_size Numeric size for tip/node points.
#' @param ylim Optional y-axis limits (unused in current implementation).
#' @param min_cells Integer minimum number of cells with VAF > 0 required to keep a variant.
#' @param tip_annot Optional data frame of tip annotations with columns `cell` and `annot`.
#' @param annot_scale Optional ggplot scale for annotation colors.
#' @param feature_mat Optional matrix of features (rows = features, columns = cells) to plot as a
#'   heatmap, or a *list* of such matrices to draw as separate stacked heatmap panels (e.g. one
#'   panel per level of a grouping). A single matrix behaves exactly as before.
#' @param feature_limits Numeric vector of length 2 giving limits for feature heatmap color scale.
#' @param feature_scale Optional ggplot scale for feature heatmap colors. When `feature_mat` is a
#'   list, either one scale (recycled across panels, so they share a colour mapping) or a list of
#'   scales, one per panel.
#' @param rescale Logical; if `TRUE`, z-score features per row before plotting.
#' @param title Optional plot title.
#' @param ytitle Optional VAF heatmap y-axis title.
#' @param xtitle Optional VAF heatmap x-axis title.
#' @param label_site Logical; included for compatibility with earlier calls.
#' @param cell_annot Optional data frame or list of data frames for annotation bars.
#' @param tip_lab Logical; if `TRUE`, show tip labels on the tree.
#' @param node_lab Logical; if `TRUE`, show internal node labels on the tree.
#' @param layered Logical; if `TRUE`, render annotation bars as layered tiles.
#' @param annot_bar_height Numeric height for each annotation bar panel.
#' @param clade_bar_height Numeric height for clade bar panel (unused in current implementation).
#' @param feature_height Numeric height for the feature heatmap panel; when `feature_mat` is a list,
#'   either one value (recycled) or one per panel.
#' @param het_max Numeric maximum VAF for heatmap color scaling.
#' @param conf_min Numeric minimum value for confidence color scale.
#' @param conf_max Numeric maximum value for confidence color scale.
#' @param conf_label Logical; if `TRUE`, label node confidence values.
#' @param branch_length Logical; if `FALSE`, drop branch lengths before plotting.
#' @param node_conf Logical; if `TRUE`, plot node confidence (requires `tbl_graph`).
#' @param annot_pal Optional palette (vector or list) for annotation bars.
#' @param annot_legend Logical; if `TRUE`, show annotation legends.
#' @param label_group Logical; if `TRUE`, label groups in annotation bars.
#' @param text_size Numeric base text size for labels.
#' @param annot_title_size Numeric size for annotation strip titles.
#' @param node_label_size Numeric size for node/tip labels.
#' @param mut Optional mutation column name in `tree` to color nodes by VAF.
#' @param mark_low_cov Logical; if `TRUE`, mark low-coverage cells on the VAF heatmap.
#' @param facet_by_group Logical; if `TRUE`, facet VAF heatmap by `group` column.
#' @param flip Logical; if `TRUE`, flip the tree orientation.
#' @param ladderize Logical; if `TRUE`, ladderize the tree.
#' @param node_scores Optional named numeric vector of node/tip scores used to color tree branches.
#'   Names may be tip labels, internal node labels, or numeric node IDs.
#' @param node_score_limits Numeric vector of length 2 giving limits for node score colors.
#' @param variants_highlight Optional vector of variant names to bold on the heatmap axis.
#' @param show_variant_names Logical; if `TRUE`, show variant tick labels on the VAF heatmap; hide when `FALSE`.
#' @param show_tree_y_axis Logical; if `TRUE`, show y-axis ticks/labels on the tree panel.
#' @param feature_legend Logical; if `TRUE`, display the legend for the feature heatmap; otherwise hide it.
#' @param raster Logical; if `TRUE`, rasterize each plot panel via `ggrastr`.
#' @param raster_dpi Numeric DPI to use when rasterizing panels.
#'
#' @return A patchwork/ggplot object combining the tree, annotations, and heatmaps.
#' @export 
plot_phylo_heatmap2 = function(phylo, df_var = NULL, branch_width = 0.25, root_edge = TRUE, dot_size = 1, ylim = NULL, min_cells = 1,
    tip_annot = NULL, annot_scale = NULL, feature_mat = NULL, feature_limits = c(-2,2), feature_scale = NULL, rescale = FALSE,
    title = NULL, ytitle = NULL, xtitle = NULL, label_site = FALSE, cell_annot = NULL, tip_lab = FALSE, node_lab = FALSE, layered = FALSE, annot_bar_height = 0.1, clade_bar_height = 1, feature_height = 1,
    het_max = 0.1, conf_min = 0, conf_max = 1, conf_label = FALSE, branch_length = TRUE, node_conf = FALSE, annot_pal = NULL, annot_legend = FALSE, label_group = FALSE,
    text_size = 3, annot_title_size = text_size, node_label_size = 1, mut = NULL, mark_low_cov = FALSE, facet_by_group = FALSE, flip = TRUE, ladderize = TRUE,
    node_scores = NULL, node_score_limits = c(-4,4), variants_highlight = NULL, show_variant_names = TRUE, show_tree_y_axis = FALSE, feature_legend = TRUE,
    raster = FALSE, raster_dpi = 300) {
    
    if (!branch_length) {
        phylo$edge.length = NULL
    }
    
    p_tree = phylo %>%
        add_node_names() %>%
        ggtree(ladderize = ladderize, linewidth = branch_width, right = flip) + 
        theme_bw() +
        theme(
            plot.margin = margin(l = 0, r = 0, t = 0, b = 0.25, unit = "mm"), 
            axis.title.x = element_blank(), 
            axis.line.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(), 
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.ticks.length.x = unit(0, "pt"),
            panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid = element_blank()
        ) + 
        coord_flip() +
        scale_x_reverse(expand = expansion(mult = 0.05)) +
        scale_y_continuous(expand = expansion(add = 1)) +
        ggtitle(title)

    if (show_tree_y_axis) {
        p_tree = p_tree + theme(axis.text.y = element_text(size = 6, hjust = 1), axis.ticks.y = element_line())
    }

    if (!is.null(node_scores)) {
        branch_df <- data.frame(node = names(node_scores), score = unname(node_scores)) %>%
            mutate(node = as.integer(node))
        p_tree <- p_tree %<+% branch_df + 
            aes(color = score) +
            scale_color_gradient2(low = "blue", mid = 'gray60', high = "red", limits = node_score_limits, oob = scales::oob_squish, na.value = 'gray60')
        p_tree <- p_tree + ggnewscale::new_scale_color()
    }

    # plot inferred ancestral VAF
    if (!is.null(mut)) {

        gtree = phylo_to_gtree(phylo)

        dat = gtree %>% activate(nodes) %>%
            as.data.frame() %>% select(all_of(c('name', p_v = mut)))

        p_tree = p_tree %<+% 
            dat + 
            ggraph::geom_node_point(aes(color = p_v), size = dot_size, stroke = 0) +
            scale_color_gradient(limits = c(0,het_max), oob = scales::oob_squish)

    }

    if (node_conf) {

        gtree = phylo_to_gtree(phylo) %>% mutate(conf = as.numeric(label))

        dat = gtree %>% activate(nodes) %>%
            mutate(isRoot = node_is_root()) %>%
            filter(!is.na(conf)) %>%
            as.data.frame() %>% 
            select(any_of(c('name', 'isRoot', 'conf')))

        if ('conf' %in% colnames(dat)) {

            p_tree = p_tree %<+% 
                dat + 
                geom_nodepoint(aes(color = conf, subset = !isTip & !isRoot, x = branch), size = dot_size, pch = 16, stroke = 0) +
                scale_color_gradient(name = 'Conf', low = 'white', high = 'firebrick', limits = c(conf_min, conf_max), oob = scales::oob_squish)

            if (conf_label) {
                p_tree = p_tree + 
                    geom_text2(
                        aes(label = round(conf, 2), x = branch, subset = !isTip & !isRoot),
                        size = text_size, hjust = 0, vjust = 0.25
                    )
            }
        }

    }

    if (!is.null(tip_annot)) {

        dat = tip_annot %>% select(any_of(c('name' = 'cell', 'annot')))

        p_tree = p_tree %<+% 
            dat +
            geom_tippoint(aes(color = annot, subset = !is.na(annot)), size = dot_size, pch = 19, stroke = 0, show.legend = FALSE)
    }

    if (tip_lab) {
        p_tree = p_tree + geom_tiplab(size = node_label_size, vjust = 0.5, hjust = 1.2, angle = 90) +
            scale_x_reverse(expand = expansion(mult = c(0.15, 0.05)))
    }

    if (node_lab) {
        p_tree = p_tree + 
            geom_text2(
                aes(label = label, subset = !isTip, x = branch),
                size = node_label_size, vjust = 0.5, hjust = 0.5
            )
    }

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% 
        pull(label)

    p_heatmap <- NULL

    if (!is.null(df_var)) {

        if (!'vaf' %in% colnames(df_var)) {
            df_var = df_var %>% mutate(vaf = a/d)
        }

        df_var = df_var %>% filter(cell %in% cell_order)

        mut_order = order_muts(cell_order, df_var)
        
        df_var = df_var %>%
            group_by(variant) %>%
            filter(sum(vaf>0, na.rm = TRUE)>=min_cells) %>%
            ungroup() %>%
            mutate(variant = factor(variant, rev(mut_order))) %>%
            mutate(cell = factor(as.integer(factor(cell, cell_order)), 1:length(cell_order)))

        xtitle = ifelse(is.null(xtitle), paste0('Cells (n=', length(cell_order), ')'), xtitle)
        ytitle = ifelse(is.null(ytitle), paste0('Variants (n=', length(unique(df_var$variant)), ')'), ytitle)

        p_heatmap = df_var %>% 
            ggplot(
                aes(x = cell, y = variant, fill = vaf)
            ) +
            geom_raster() +
            theme_bw() + 
            theme(
                plot.margin = margin(t = 0.25, r = 0, b = 0, l = 0, unit = "mm"), 
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                panel.background = element_rect(fill = "transparent", colour = NA), 
                plot.background = element_rect(fill = "transparent"),
                panel.grid = element_blank()
            ) +
            scale_x_discrete(expand = expansion(add = 1), drop = F) +
            scale_y_discrete(expand = expansion(mult = 0.01)) +
            scale_fill_gradient(low = 'white', high = 'red', limits = c(0,het_max), oob = scales::oob_squish) +
            guides(fill = guide_colorbar(title = 'VAF')) +
            xlab(xtitle) +
            ylab(ytitle)

        if (show_variant_names) {
            p_heatmap = p_heatmap + theme(axis.text.y = element_text(size = text_size))
        } else {
            p_heatmap = p_heatmap + theme(axis.text.y = element_blank())
        }

        if (!is.null(variants_highlight)) {
            vh <- intersect(as.character(unique(df_var$variant)), variants_highlight)
            if (length(vh) > 0) {
                lab_expr <- function(v) sapply(v, function(x) if (x %in% vh) paste0('bold("', x, '")') else paste0('"', x, '"'))
                p_heatmap <- p_heatmap +
                    scale_y_discrete(labels = function(v) parse(text = lab_expr(v)))
            }
        }

        if (facet_by_group) {
            p_heatmap = p_heatmap + facet_grid(group~., scales = 'free_y', space = 'free_y') +
                theme(panel.spacing.y = unit(0,'mm'))
        }
        
        if (mark_low_cov) {
            p_heatmap = p_heatmap + geom_point(
                data = df_var %>% filter(d < 10),
                pch = 4,
                size = dot_size
            )
        }

        if (raster) {
            p_heatmap = ggrastr::rasterize(p_heatmap, dpi = raster_dpi, layers = 'Raster')
        }
    }

    if (!is.null(cell_annot)) {
        
        # Handle both single data frame and list of data frames
        if (!is.list(cell_annot) || is.data.frame(cell_annot)) {
            cell_annot <- list(cell_annot)
        }

        # Handle annot_pal - if single value, replicate for all annotations
        if (!is.null(annot_pal)) {
            if (!is.list(annot_pal)) {
                annot_pal <- rep(list(annot_pal), length(cell_annot))
            }
        } else {
            annot_pal <- rep(list(NULL), length(cell_annot))
        }

        # Filter each annotation to only include cells in the tree
        cell_annot <- lapply(cell_annot, function(annot) {
            annot %>% filter(cell %in% phylo$tip.label)
        })

        if (is.null(names(cell_annot))) {
            names(cell_annot) <- paste0('annot', seq_along(cell_annot))
        }
        annot_titles <- names(cell_annot)

        # Create separate bars for each annotation with its own palette
        p_bars <- mapply(function(annot, pal, bar_title, legend_title) {
            annot %>%
                mutate(cell = factor(cell, cell_order)) %>%  
                annot_bar(
                    legend = annot_legend, 
                    label_group = label_group, 
                    label_size = annot_title_size,
                    annot_pal = pal,
                    annot_scale = annot_scale,
                    legend_title = legend_title, 
                    layered = layered,
                    raster = raster,
                    raster_dpi = raster_dpi
                ) +
                theme(
                    plot.margin = margin(t = 0.25, r = 0, b = 0.25, l = 0, unit = 'mm'),
                    panel.border = element_rect(size = 0.25, color = 'black', fill = NA)
                )
        }, cell_annot, annot_pal, annot_titles, annot_titles, SIMPLIFY = FALSE)

    }

    if (!is.null(feature_mat)) {

        # Accept either a single matrix or a list of matrices (drawn as separate stacked panels).
        # Mirrors how `cell_annot` is handled above. A bare matrix is not a list, but a data.frame
        # is, hence the second test.
        if (!is.list(feature_mat) || is.data.frame(feature_mat)) {
            feature_mat <- list(feature_mat)
        }
        n_feat <- length(feature_mat)

        if (is.null(feature_scale)) {
            feature_scale_list <- rep(list(NULL), n_feat)
        } else if (inherits(feature_scale, "Scale") || !is.list(feature_scale)) {
            feature_scale_list <- rep(list(feature_scale), n_feat)   # one scale shared by all panels
        } else {
            feature_scale_list <- feature_scale
        }
        feature_height_vec <- rep_len(feature_height, n_feat)

        p_feature_list <- lapply(seq_len(n_feat), function(fi) {

            fm <- feature_mat[[fi]]

            df_feature = fm %>%
                as.data.frame() %>%
                tibble::rownames_to_column('feature') %>%
                mutate(feature = factor(feature, rev(rownames(fm)))) %>%
                reshape2::melt(id.vars = 'feature', variable.name = 'cell', value.name = 'value') %>%
                mutate(cell = factor(cell, cell_order))

            if (rescale) {
                df_feature = df_feature %>%
                    group_by(feature) %>%
                    mutate(value = as.vector(scale(value))) %>%
                    ungroup()
            }

            p_feature = df_feature %>%
                ggplot(aes(x = cell, y = feature, fill = value)) +
                geom_raster(show.legend = feature_legend) +
                theme_bw() +
                theme(axis.text.x = element_blank(),
                    axis.text.y = element_text(size = text_size),
                    plot.margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "mm"),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                ) +
                scale_x_discrete(expand = expansion(add = 1), drop = F) +
                scale_y_discrete(expand = expansion(add = 0))

            if (is.null(feature_scale_list[[fi]])) {
                p_feature = p_feature + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits = feature_limits, oob = scales::oob_squish)
            } else {
                p_feature = p_feature + feature_scale_list[[fi]]
            }

            if (raster) {
                p_feature = ggrastr::rasterize(p_feature, dpi = raster_dpi, layers = 'Raster')
            }
            p_feature
        })
    }

    # Determine heights based on components
    heights <- c(1)  # tree always has height 1

    # Build plot components list
    plot_components <- list(p_tree)
    
    if (!is.null(cell_annot)) {
        plot_components <- c(plot_components, p_bars)
        heights <- c(heights, rep(annot_bar_height, length(p_bars)))
    }
    
    if (!is.null(feature_mat)) {
        plot_components <- c(plot_components, p_feature_list)
        heights <- c(heights, feature_height_vec)
    }
    
    if (!is.null(p_heatmap)) {
        plot_components <- c(plot_components, list(p_heatmap))
        heights <- c(heights, 2)
    }

    # Combine plots
    wrap_plots(plot_components) + plot_layout(heights = heights, guides = 'collect')
}

#' Create an annotation bar plot
#'
#' Builds a ggplot annotation bar (horizontal colored tiles) from a data
#' frame with `cell` and `annot` columns.
#'
#' @param D Data frame with at least columns `cell` and `annot`.
#' @param transpose Logical; if `TRUE`, flip the bar orientation.
#' @param legend Logical; if `TRUE`, show fill legend.
#' @param legend_title Character legend title.
#' @param size Numeric; not currently used (retained for API compatibility).
#' @param label_group Logical; if `TRUE`, add text labels for each group.
#' @param label_size Numeric text size for group labels and axis text.
#' @param annot_pal Optional named color vector for `scale_fill_manual`.
#' @param annot_scale Optional ggplot fill scale override.
#' @param raster Logical; if `TRUE`, rasterize via `ggrastr`.
#' @param raster_dpi Numeric DPI for rasterization.
#' @param layered Logical; if `TRUE`, plot annotations as layered tiles
#'   instead of a single row.
#' @return A ggplot object.
#' @keywords internal
#' @noRd
annot_bar = function(
    D, transpose = FALSE, legend = TRUE, legend_title = '', size = 0.05, label_group = FALSE, label_size = 5,
    annot_pal = NULL, annot_scale = NULL, raster = FALSE, raster_dpi = 300, layered = FALSE
) {

    if (layered) {
        p = ggplot(D, aes(x = cell, y = annot, fill = annot)) +
            scale_y_discrete(expand = expansion(add = 1.5))
    } else {
        p = ggplot(D, aes(x = cell, y = legend_title, fill = annot)) +
            scale_y_discrete(expand = expansion(mult = 0.01))
    }

    p = p +
        geom_raster() +
        theme_void() +
        scale_x_discrete(expand = expansion(add = 1), drop = F) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0, color = 'black', fill = NA),
            panel.background = element_rect(fill = '#CACACA'),
            strip.background = element_blank(),
            strip.text = element_blank(),
            # axis.text = element_text(size = 8),
            axis.text.y = element_text(size = label_size, hjust = 1, vjust = 0.5, margin = margin(r = 0.5, unit = 'mm')),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(0.1,0,0.1,0, unit = 'mm')
        )

    if (label_group) {
        p = p + geom_text(aes(label = annot), size = label_size, angle = 90)
    }


    if (!is.null(annot_pal)) {
        p = p + scale_fill_manual(values = annot_pal, na.value = 'gray90', limits = force)
    } else {
        p = p + scale_fill_discrete(na.value = 'gray90', limits = force)
    }

    if (!is.null(annot_scale)) {
        p = p + annot_scale
    }

    if (transpose) {
        p = p + coord_flip() +
            theme(plot.margin = margin(0,0.5,0,0.5, unit = 'mm'))
    }

    if (legend) {
        p = p + guides(fill = guide_legend(keywidth = unit(3, 'mm'), keyheight = unit(1, 'mm'), title = legend_title))
    } else {
        p = p + guides(fill = 'none')
    }

    if (raster) {
        p = ggrastr::rasterize(p, layers = 'Raster', dpi = raster_dpi)
    }

    return(p)
}

#' Order mutations by weighted-median tip position
#'
#' Sorts variants by the VAF-weighted median tip position, producing a
#' heatmap ordering that aligns with the tree layout.
#'
#' @param cell_order Character vector of cell names in tree-tip order.
#' @param mut_dat Data frame with columns `variant`, `cell`, and `vaf`.
#' @return Character vector of variant names in sorted order.
#' @keywords internal
#' @noRd
order_muts <- function(cell_order, mut_dat) {
  
  # 1) cast to a variant × cell matrix of continuous VAFs
  vaf_mat <- mut_dat %>% 
    reshape2::dcast(variant ~ cell, value.var = 'vaf') %>%
    tibble::column_to_rownames('variant')
  
  # 2) make sure columns follow the tree’s tip order
  vaf_mat <- vaf_mat[, cell_order, drop = FALSE]
  
  # 3) create an index for each tip
  tip_idx <- seq_along(cell_order)
  names(tip_idx) <- cell_order
  
  # 4) compute weighted-median tip for each variant
  med_tip_pos <- apply(vaf_mat, 1, function(vafs) {
    # treat missing as zero
    vafs[is.na(vafs)] <- 0
    tot <- sum(vafs)
    if (tot == 0) {
      return(NA_real_)    # no signal → put at end
    }
    # cumulative weight
    cs <- cumsum(vafs)
    # first tip where cum-weight ≥ half total
    k  <- which(cs >= tot/2)[1]
    tip_idx[k]
  })
  
  # 5) sort variants by that median tip
  names(sort(med_tip_pos, na.last = TRUE))
}

#' Plot Circular Phylogenetic Tree with Annotations
#'
#' Creates a circular (fan) layout phylogenetic tree plot with optional annotations,
#' feature heatmaps, and node confidence values.
#'
#' @param gtree A `phylo` object or `tbl_graph` object representing the tree.
#' @param cell_annot A data frame or list of data frames containing cell annotations.
#'   Columns should include 'cell' and 'annot'.
#' @param tip_annot A data frame containing tip annotations. Columns should include 'cell' and 'annot'.
#' @param feature_mat A matrix of features to plot as a heatmap (rows are features, columns are cells).
#' @param branch_scores Optional branch scores (currently unused).
#' @param node_conf Logical. If `TRUE` and `gtree` is a `tbl_graph`, plots node confidence.
#' @param conf_label Logical. If `TRUE`, adds text labels for node confidence.
#' @param title Character string for the plot title.
#' @param pwidth_annot Numeric. Width of the annotation ring(s).
#' @param pwidth_feature Numeric. Width of the feature heatmap ring.
#' @param branch_width Numeric. Line width for tree branches.
#' @param dot_size Numeric. Size of tip/node points.
#' @param conf_min Numeric. Minimum value for confidence color scale.
#' @param conf_max Numeric. Maximum value for confidence color scale.
#' @param annot_pal A vector or list of vectors specifying colors for annotations.
#' @param offset Numeric. Offset distance between tree and annotation rings.
#' @param width Numeric. Width of heatmap cells.
#' @param label_size Numeric. Size of text labels.
#' @param rescale Logical. If `TRUE`, rescales features (z-score) for the heatmap.
#' @param limits Numeric vector of length 2. Limits for the feature heatmap color scale.
#' @param feature_legend_title Character string. Title for the feature heatmap legend.
#' @param flip Logical. If `TRUE`, flips the tree direction.
#' @param legend Logical. If `TRUE`, shows legends.
#' @param layered Logical. If `TRUE`, plots annotations as layered tiles.
#' @param smooth_k Integer. Number of neighbors for smoothing features. 0 means no smoothing.
#' @param ladderize Logical. If `TRUE`, ladderizes the tree.
#' @param open_angle Numeric. Angle of the opening in the fan layout.
#' @param feature_axis_label Logical. If `TRUE`, show feature heatmap axis labels; hide when `FALSE`.
#' @param feature_axis_angle Numeric. Rotation angle (degrees) for feature axis text.
#'
#' @return A `ggtree` plot object.
#' @importFrom ggtreeExtra geom_fruit
#' @export
plot_phylo_circ = function(gtree, cell_annot = NULL, tip_annot = NULL, feature_mat = NULL, branch_scores = NULL,
    node_conf = FALSE, conf_label = FALSE, title = '', pwidth_annot = 0.25, pwidth_feature = 0.25,
    branch_width = 0.3, dot_size = 1, conf_min = 0, conf_max = 0.5, annot_pal = NULL, offset = 0.05, width = 0.8,
    label_size = 2, rescale = FALSE, limits = c(-2,2), feature_legend_title = 'Score', flip = TRUE,
    annot_legend = TRUE, feature_legend = TRUE, layered = FALSE, smooth_k = 0, ladderize = TRUE, open_angle = 0, feature_axis_label = TRUE,
    feature_axis_angle = 30) {

    p_tree = ggtree(gtree, 
            ladderize = ladderize, layout = 'fan', open.angle = open_angle,
            branch.length = "none", linewidth = branch_width, right = flip
        ) + ggtitle(title)

    if (node_conf & inherits(gtree, 'tbl_graph')) {

        dat = gtree %>% activate(nodes) %>%
            mutate(isRoot = node_is_root()) %>%
            as.data.frame() %>% 
            select(any_of(c('name', 'isRoot', 'conf')))

        if ('conf' %in% colnames(dat)) {
            p_tree = p_tree %<+% 
                dat +
                geom_nodepoint(aes(color = conf, subset = !isTip & !isRoot, x = branch), size = dot_size, pch = 16, stroke = 1) +
                scale_color_gradient(low = 'white', high = 'firebrick', limits = c(conf_min, conf_max), oob = scales::oob_squish)
    
            if (conf_label) {
                p_tree = p_tree + 
                    geom_text2(
                        aes(label = round(conf, 2), x = branch, subset = !isTip & !isRoot),
                        size = text_size, hjust = 0, vjust = 0.25
                    )
            }
        }
    }

    if (!is.null(tip_annot)) {

        dat = tip_annot %>% select(any_of(c('name' = 'cell', 'annot')))

        p_tree = p_tree %<+% 
            dat +
            geom_tippoint(aes(color = annot), size = dot_size, pch = 19, stroke = 0)
    }

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% 
        pull(label)


    if (!is.null(cell_annot)) {
        
        # Handle both single data frame and list of data frames
        if (!is.list(cell_annot) || is.data.frame(cell_annot)) {
            cell_annot <- list(cell_annot)
        }
        
        # Handle annot_pal - if single value, replicate for all annotations
        if (!is.null(annot_pal)) {
            if (!is.list(annot_pal)) {
                annot_pal <- rep(list(annot_pal), length(cell_annot))
            }
        } else {
            annot_pal <- rep(list(NULL), length(cell_annot))
        }
        
        # Ensure annotation list has names to use as legend titles
        if (is.null(names(cell_annot))) {
            names(cell_annot) <- paste0('annot', seq_along(cell_annot))
        }

        # Add each annotation as a separate fruit layer
        for (i in seq_along(cell_annot)) {
            annot_data <- cell_annot[[i]]
            pal <- annot_pal[[i]]
            bar_title <- names(cell_annot)[i]
            
            # Start new scale for each annotation
            p_tree = p_tree + ggnewscale::new_scale_fill()
            
            if (layered) {
                p_tree = p_tree + geom_fruit(
                    data = annot_data,
                    geom = geom_tile,
                    pwidth = pwidth_annot,
                    offset = offset,
                    mapping = aes(y = cell, fill = annot, x = annot),
                    show.legend = annot_legend)
            } else {
                p_tree = p_tree + geom_fruit(
                    data = annot_data,
                    geom = geom_col,
                    pwidth = pwidth_annot,
                    offset = offset,
                    mapping = aes(y = cell, fill = annot, x = 1),
                    show.legend = annot_legend)
            }
            
            # Add custom palette if provided; use the per-annotation name as legend title
            na_col = 'gray60'
            if (!is.null(pal)) {
                p_tree = p_tree + scale_fill_manual(name = bar_title, values = pal, na.value = na_col)
            } else {
                p_tree = p_tree + scale_fill_discrete(name = bar_title, na.value = na_col)
            }
        }
    }

    if (!is.null(feature_mat)) {

        feature_mat = as.matrix(feature_mat)

        feature_mat = feature_mat[,cell_order,drop = FALSE] 

        # if (rescale) {
        #     feature_mat = t(scale(t(feature_mat)))
        # }

        if (smooth_k > 0) {
            feature_mat = feature_mat %>% row_smooth(k = smooth_k)
        }

        df_feature = feature_mat %>%
            reshape2::melt() %>% 
            setNames(c('feature', 'cell', 'value')) %>%
            mutate(cell = as.character(cell))

        if (rescale) {
            df_feature = df_feature %>%
                group_by(feature) %>%
                mutate(value = as.vector(scale(value))) %>%
                ungroup()
        }

        axis_flag = if (feature_axis_label) "x" else "none"
        axis_text_size = if (feature_axis_label) label_size else 0

        p_tree = p_tree + 
            ggnewscale::new_scale_fill() +
            geom_fruit(
                data = df_feature,
                geom = geom_tile,
                offset = offset,
                pwidth = pwidth_feature,
                width = width,
                # linewidth = 0.3,
                axis.params = list(
                    axis       = axis_flag, 
                    text.angle = feature_axis_angle,
                    text.size  = axis_text_size,
                    vjust      = 0.5,
                    hjust = 1
                ),
                mapping = aes(x = feature, y = cell, fill = value),
                show.legend = feature_legend
            ) + 
            scale_fill_gradient2(name = feature_legend_title, low = "blue", high = "red", limits = limits, oob = scales::oob_squish)

        # create fake legend for tip annotation
        # if (!is.null(tip_annot)) {
        #     cts = unique(tip_annot$annot)
        #     n_cts = length(cts)
        #     max_set1 = RColorBrewer::brewer.pal.info['Set1','maxcolors']
        #     if (n_cts <= max_set1) {
        #         annot_cols = RColorBrewer::brewer.pal(n_cts, 'Set1') %>% setNames(cts)
        #     } else {
        #         base_cols = RColorBrewer::brewer.pal(max_set1, 'Set1')
        #         annot_cols = grDevices::colorRampPalette(base_cols)(n_cts) %>% setNames(cts)
        #     }

        #     p_tree = p_tree +
        #         scale_color_manual(
        #             name   = "Annotation",   
        #             values = annot_cols
        #         ) +
        #         geom_point(
        #             data       = data.frame(annot = names(annot_cols)),
        #             aes(colour  = annot),
        #             x          = 0,      # plotted off the panel
        #             y          = 0,
        #             size = 0, stroke = 0,
        #             show.legend= TRUE
        #         )  +
        #         guides(color = guide_legend(override.aes = list('size' = 2), ncol = 2)) 
        # }
    }

    return(p_tree)
}


##
#' Moving-average smoothing across matrix rows
#'
#' Applies a fixed-width moving average independently to every row of a matrix.
#' Missing values can be optionally ignored, and the `edge` argument controls how
#' windows at the matrix boundaries are handled. With `edge = "partial"`
#' (default), windows that extend past the matrix edge shrink to the available
#' cells, while `edge = "full"` enforces the full window size by returning
#' `NA` when a complete window cannot be formed.
#'
#' @param M Numeric matrix to smooth (rows are independent series).
#' @param k Integer window width (must be >= 1).
#' @param na_rm Logical; if `TRUE`, compute means ignoring `NA`s, otherwise only
#'   windows containing no `NA` values contribute.
#' @param edge Either ``"partial"`` or ``"full"``; governs boundary handling as
#'   described above.
#'
#' @return Matrix of the same dimensions as `M` containing smoothed values.
#' @keywords internal
row_smooth <- function(M, k, na_rm = FALSE, edge = c("partial", "full")) {
	edge <- match.arg(edge)
	if (!is.matrix(M)) stop("M must be a matrix.")
	if (k < 1L) stop("k must be >= 1.")
	if (k == 1L) return(M)

	nr <- nrow(M)
	nc <- ncol(M)
	out <- matrix(NA_real_, nr, nc, dimnames = dimnames(M))

	left <- k %/% 2L
	right <- k - left - 1L

	idx <- seq_len(nc)
	L <- pmax(1L, idx - left)
	R <- pmin(nc, idx + right)
	win_len <- R - L + 1L

	for (i in seq_len(nr)) {
		x <- M[i, ]
		val <- ifelse(is.na(x), 0, x)
		cs <- c(0, cumsum(val))
		cc <- c(0, cumsum(!is.na(x)))

		sums <- cs[R + 1L] - cs[L]
		cnt  <- cc[R + 1L] - cc[L]

		if (na_rm) {
			den <- cnt
			if (edge == "full") den[win_len < k] <- NA_integer_
			res <- sums / den
			res[!is.finite(res)] <- NA_real_
		} else {
			ok <- (cnt == win_len)
			if (edge == "full") ok <- ok & (win_len == k)
			res <- sums / win_len
			res[!ok] <- NA_real_
		}
		out[i, ] <- res
	}
	out
}


#' Build a named Brewer palette with optional cycling behavior
#'
#' Creates up to `n` visually distinct colors by sampling (and, when needed,
#' interpolating) from an RColorBrewer palette, optionally reordering each
#' block in a zig-zag pattern and progressively lightening subsequent cycles.
#'
#' @param n Integer number of colors to generate.
#' @param labels Optional character vector of names to assign to the colors. If
#'   `NULL`, sequential integers are used.
#' @param pal RColorBrewer palette name provided to `RColorBrewer::brewer.pal`.
#' @param alternating Logical; when `TRUE`, reorders each chunk via a
#'   zig-zag pattern to maximize separation of adjacent hues.
#' @param cycle_len Positive integer forcing palette generation in repeating
#'   chunks of this size; `0` (default) disables chunking.
#' @param cycle_shift Scalar in `[0, 1]` controlling how strongly later chunks
#'   blend toward white to maintain distinguishability.
#'
#' @return Named character vector of hex colors.
#' @export
make_clade_pal <- function(
	n,
	labels = NULL,
	pal = "Set3",
	alternating = FALSE,
	cycle_len = 0L,
	cycle_shift = 0.15
) {

	max_n <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
	base_cols <- RColorBrewer::brewer.pal(max_n, pal)

	zigzag_idx <- function(len) {
		if (len <= 2L) return(seq_len(len))
		ord <- as.vector(rbind(seq_len(len), rev(seq_len(len))))
		ord <- ord[ord <= len]
		ord[!duplicated(ord)]
	}

	blend_with_white <- function(cols, frac) {
		if (frac <= 0) return(cols)
		frac <- pmin(frac, 0.9)
		rgb_mat <- grDevices::col2rgb(cols) / 255
		mix <- (1 - frac) * rgb_mat + frac
		grDevices::rgb(mix[1, ], mix[2, ], mix[3, ])
	}

	if (cycle_len <= 0L) {
		cols <- if (n <= max_n) {
			base_cols[seq_len(n)]
		} else {
			grDevices::colorRampPalette(base_cols)(n)
		}
		if (alternating) {
			cols <- cols[zigzag_idx(length(cols))]
		}
	} else {
		chunk_size <- max(1L, min(cycle_len, max_n))
		base_chunk <- if (chunk_size <= max_n) {
			base_cols[seq_len(chunk_size)]
		} else {
			grDevices::colorRampPalette(base_cols)(chunk_size)
		}

		cols <- character(0)
		remaining <- n
		cycle_idx <- 1L
		while (remaining > 0L) {
			take <- min(chunk_size, remaining)
			block <- base_chunk
			if (take < length(block)) {
				block <- block[seq_len(take)]
			}
			if (alternating && take > 1L) {
				block <- block[zigzag_idx(take)]
			}
			if (cycle_idx > 1L) {
				frac <- min(0.85, (cycle_idx - 1L) * cycle_shift)
				block <- blend_with_white(block, frac)
			}
			cols <- c(cols, block)
			remaining <- remaining - take
			cycle_idx <- cycle_idx + 1L
		}
	}

	if (is.null(labels)) {
		labels <- as.character(seq_len(n))
	} else {
		labels <- as.character(labels)
	}

	names(cols) <- labels
	cols
}
