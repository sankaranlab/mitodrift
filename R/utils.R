
#' Convert a long format mutation data frame to a matrix
#' @param mut_dat_long A data frame with columns: variant, cell, variable
#' @param variable The variable to convert to a matrix
#' @return A matrix of values
#' @export
long_to_mat = function(mut_dat_long, variable) {
    data.table::as.data.table(mut_dat_long) %>%
        data.table::dcast(variant ~ cell, value.var = variable, fill = 0) %>%
        tibble::column_to_rownames("variant") %>%
        as.matrix()
}

#' Convert allele count matrix (amat) and total count matrix (dmat) to long format
#' @param amat Allele count matrix with variants as rows and cells as columns
#' @param dmat Total count matrix with variants as rows and cells as columns
#' @return A data.table in long format with columns: variant, cell, a (allele count), d (total count)
#' @export
mat_to_long = function(amat, dmat) {
    # Ensure both matrices have the same dimensions and row/column names
    if (!identical(dim(amat), dim(dmat)) || 
        !identical(rownames(amat), rownames(dmat)) || 
        !identical(colnames(amat), colnames(dmat))) {
        stop("amat and dmat must have identical dimensions and row/column names")
    }
    
    # Convert to data.table and melt to long format
    amat_dt <- data.table::as.data.table(amat, keep.rownames = "variant")
    dmat_dt <- data.table::as.data.table(dmat, keep.rownames = "variant")
    
    # Melt both matrices to long format
    amat_long <- data.table::melt(amat_dt, id.vars = "variant", 
                                  variable.name = "cell", value.name = "a")
    dmat_long <- data.table::melt(dmat_dt, id.vars = "variant", 
                                  variable.name = "cell", value.name = "d")
    
    # Merge the two long tables
    result <- merge(amat_long, dmat_long, by = c("variant", "cell"))
    
    # Set column order
    data.table::setcolorder(result, c("variant", "cell", "a", "d"))
    
    return(result)
}

#' Find the root node ID of a rooted phylogeny.
#'
#' The root is identified as the node index that appears as a parent in
#' `tree$edge[, 1]` but never as a child in `tree$edge[, 2]`.
#'
#' @param tree A rooted `phylo` object.
#' @return Integer node ID for the root.
#' @keywords internal
#' @noRd
find_root = function(tree) {
	if (!ape::is.rooted(tree)) {
		stop("tree must be rooted")
	}
	setdiff(tree$edge[, 1], tree$edge[, 2])[1]
}

#' Get the clade node IDs adjacent to the root.
#'
#' Returns the child node IDs of the root (the clades immediately below the
#' root). For a rooted binary tree this is typically length 2.
#'
#' @param tree A rooted `phylo` object.
#' @return Integer vector of root child node IDs.
#' @keywords internal
#' @noRd
find_root_split = function(tree) {
	root = find_root(tree)
	root_children = tree$edge[tree$edge[, 1] == root, 2]
	return(root_children)
}

#' Collapse weak clades below a confidence threshold.
#'
#' Given a fully binary phylogeny with per-node confidence scores stored in
#' `tree$node.label`, this helper collapses every internal node whose confidence
#' falls below `conf`, effectively pruning low-support clades while retaining
#' higher-confidence structure.
#'
#' @param tree A binary `phylo` object whose `node.label` vector stores
#'   posterior/confidence values for internal nodes.
#' @param conf Numeric threshold in $[0,1]$. Internal nodes with confidence below
#'   this value are collapsed.
#' @param collapse_trivial Logical; if `TRUE`, collapses the trivial
#'   root-adjacent singleton split (1 tip vs the rest) into a root polytomy.
#' @return A renumbered `phylo` object with low-confidence nodes collapsed.
#' @export
trim_tree = function(tree, conf, collapse_trivial = TRUE) {
    if (!ape::is.binary(tree)) {stop("Tree must be binary")}
    n_tip = length(tree$tip.label)
    node_confs = as.numeric(tree$node.label)
    # treat NA as no confidence
    node_confs[is.na(node_confs)] = 0
    collapse_list = which(node_confs < conf) + n_tip

	# Test if trivial split (1 tip vs rest)
    if (collapse_trivial && conf > 0) {
        root_children = find_root_split(tree)
        is_tip = root_children <= n_tip
		if (any(is_tip)) {
            internal_root_child = root_children[!is_tip]
            if (!internal_root_child %in% collapse_list) {
                collapse_list = c(collapse_list, internal_root_child)
            }
		}
	}

    tree = TreeTools::CollapseNode(tree, collapse_list)
    tree = TreeTools::Renumber(tree)

    return(tree)
}

#' Collapse branches with high expected mis-assignments.
#'
#' Computes the expected number of incorrectly assigned tips for each internal
#' branch (based on (1 - p) * clade size, where p is the branch
#' confidence) and collapses branches whose expectation exceeds a tolerance
#' derived from the total number of tips.
#'
#' @param tree A binary `phylo` object with internal node confidence scores in
#'   `tree$node.label`.
#' @param tol Numeric tolerance expressed as a fraction of total tips; the
#'   threshold is `length(tree$tip.label) * tol` expected errors.
#' @return A renumbered `phylo` object with overly uncertain branches collapsed.
#' @export
trim_tree_exp = function(tree, tol) {
	
	n_tip <- length(tree$tip.label)
    n_tol = n_tip * tol
	
	# branch posteriors on internal nodes
	node_confs <- as.numeric(tree$node.label)
	node_confs[is.na(node_confs)] <- 0
	
	internal_nodes <- (n_tip + 1L):(n_tip + tree$Nnode)
	desc_nodes <- TreeTools::CladeSizes(tree, nodes = internal_nodes)
	clade_tips <- as.integer((desc_nodes + 2L) / 2L)
	
	# Expected wrong cells per branch: (1 - p) * clade_size
	exp_wrong <- (1 - node_confs) * clade_tips

	# Collapse branches whose expected wrong cells exceeds epsilon (= n_tol)
	nodes_to_collapse <- which(exp_wrong > n_tol) + n_tip
	
	if (length(nodes_to_collapse) > 0L) {
		tree <- TreeTools::CollapseNode(tree, nodes_to_collapse)
		tree <- TreeTools::Renumber(tree)
	}
	
	return(tree)
}

#' Assign clone IDs to a tree allowing small polytomies.
#'
#' @param tree Rooted `phylo` object with bifurcating or polytomous structure.
#' @param k Positive scalar; maximum clade size that collapses into one clone.
#' @param paraphyletic Logical; if `TRUE`, leftover unassigned tips below a node
#'   are grouped into one clone, otherwise each becomes a singleton.
#' @param return_df Logical; if `TRUE`, returns a data frame annotated per tip,
#'   otherwise an integer vector of clone assignments.
#'
#' @return Either a data frame with columns `cell`, `clade`, `annot`, `size`,
#'   `frac` or an integer vector of clone IDs indexed by tip order.
#'
#' @details Finds the root (node never used as a child), recurses over internal
#' nodes, and assigns clone IDs whenever a clade has `<= k` tips. Tips directly
#' attached to the root are always singleton clones. Ensures all tips receive
#' an assignment, warning if late singletons are required.
#' @export
assign_clones_polytomy <- function(tree, k = Inf, paraphyletic = FALSE, return_df = TRUE) {

	Ntip <- length(tree$tip.label)
	if (Ntip == 0L) {
		return(integer(0L))
	}

	# Find root as the node that never appears as a child
	all_nodes <- unique(as.vector(tree$edge))
	child_nodes <- unique(tree$edge[, 2])
	root_node <- setdiff(all_nodes, child_nodes)
	if (length(root_node) != 1L) {
		stop("Could not uniquely identify root node")
	}
	root_node <- root_node[1L]

	total_nodes <- Ntip + tree$Nnode

	# children_list: parent -> vector of children node indices
	children_list <- split(tree$edge[, 2], tree$edge[, 1])

	# Precompute descendant tips for each node (memoized)
	subtree_tips <- vector("list", total_nodes)
	get_subtree_tips <- function(nd) {
		if (nd <= Ntip) return(nd)
		if (!is.null(subtree_tips[[nd]])) return(subtree_tips[[nd]])

		ch <- children_list[[as.character(nd)]]
		if (is.null(ch) || length(ch) == 0L) {
			subtree_tips[[nd]] <<- integer(0L)
			return(integer(0L))
		}
		tips_nd <- integer(0L)
		for (cnd in ch) {
			if (cnd <= Ntip) {
				tips_nd <- c(tips_nd, cnd)
			} else {
				tips_nd <- c(tips_nd, get_subtree_tips(cnd))
			}
		}
		tips_nd <- unique(tips_nd)
		subtree_tips[[nd]] <<- tips_nd
		tips_nd
	}

	# clone_assign[i] = clone ID for tip i (0 = unassigned)
	clone_assign <- integer(Ntip)
	clone_node_assign <- integer(Ntip)  # track which node defines each clone
	next_clone <- 1L

	assign_tips <- function(tips, node_id) {
		tips <- as.integer(tips)
		if (length(tips) == 0L) return(invisible(NULL))
		# only assign still-unassigned tips
		tips <- tips[clone_assign[tips] == 0L]
		if (length(tips) == 0L) return(invisible(NULL))
		clone_assign[tips] <<- next_clone
		clone_node_assign[tips] <<- node_id
		next_clone <<- next_clone + 1L
		invisible(NULL)
	}

	# Recursive processing for internal nodes below the root
	process_node <- function(nd) {
		# nd is guaranteed > Ntip here
		node_tips <- get_subtree_tips(nd)

		# If this whole clade is small enough, make it a single clone
		if (length(node_tips) <= k) {
			assign_tips(node_tips, nd)
			return(invisible(NULL))
		}

		# Otherwise, recurse into internal children first
		ch <- children_list[[as.character(nd)]]
		if (!is.null(ch) && length(ch) > 0L) {
			for (cnd in ch) {
				if (cnd > Ntip) {
					process_node(cnd)
				}
			}
		}

		# After children have formed their clones,
		# group any leftover tips under this node
		remaining <- node_tips[clone_assign[node_tips] == 0L]
		if (length(remaining) > 0L) {
			if (paraphyletic) {
				# group all remaining paraphyletic tips into a single clone
				assign_tips(remaining, nd)
			} else {
				# assign each remaining tip as its own singleton clone
				for (tip in remaining) {
					assign_tips(tip, tip)
				}
			}
		}

		invisible(NULL)
	}

	# --- 1) Handle root-level children ---
	root_children <- children_list[[as.character(root_node)]]

	# a) Tips directly from root are always their own singleton clones
	for (child in root_children) {
		if (child <= Ntip) {
			assign_tips(child, child)
		}
	}

	# b) Internal children are processed with the recursive logic
	for (child in root_children) {
		if (child > Ntip) {
			process_node(child)
		}
	}

	# Sanity: everything should be assigned
	if (any(clone_assign == 0L)) {
		warning("Some tips remain unassigned; they will be given unique singleton clones.")
		unassigned <- which(clone_assign == 0L)
		for (i in unassigned) {
			assign_tips(i, i)
		}
	}

	names(clone_assign) <- tree$tip.label
	names(clone_node_assign) <- tree$tip.label

	if (return_df) {
		clade_annot = data.frame(clade = clone_assign, clade_node = clone_node_assign) %>%
			tibble::rownames_to_column("cell") %>%
			mutate(annot = as.character(clade)) %>%
			group_by(clade) %>%
			mutate(size = n()) %>%
			ungroup() %>%
			mutate(frac = size / n()) %>%
			arrange(-frac) %>%
            mutate(clade = as.integer(factor(clade, unique(clade)))) %>%
            mutate(clade = as.character(clade)) %>%
            mutate(annot = clade)
		return(clade_annot)
	} else {
		return(list(
			clade = clone_assign,
			clade_node = clone_node_assign
		))
	}
}


#' Add "Node<n>" labels to the internal nodes of a phylo tree
#'
#' @param tree A phylo object
#' @param prefix Character prefix for node names (default "Node")
#' @param start_from_tip Logical; if `TRUE` (default), node numbering starts
#'   from `ntip + 1` (standard ape convention); if `FALSE`, numbering starts
#'   from 1.
#' @return The same phylo object, with tree$node.label set to prefix + node numbers
#' @export
add_node_names <- function(tree, prefix = "Node", start_from_tip = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylo object")
  }
  ntip  <- length(tree$tip.label)
  nnode <- tree$Nnode

  # internal node IDs run from ntip+1 to ntip+nnode
  if (start_from_tip) {
    ids <- seq(ntip + 1, ntip + nnode)
  } else {
    ids <- seq(1, nnode)
  }

  # assign labels
  tree$node.label <- paste0(prefix, ids)

  return(tree)
}


#' Convert a `phylo` object to a `tbl_graph`
#'
#' Builds a `tbl_graph` with tip nodes labeled by `phy$tip.label` and internal
#' nodes labeled from `phy$node.label`. If internal labels are missing, they are
#' created via `add_node_names()`.
#'
#' @param phy A `phylo` object.
#' @return A `tbl_graph` with `nodes` and `edges` from the input tree.
#' @export
phylo_to_gtree = function(phy) {
        
    tip_nodes <- data.frame(
        name = phy$tip.label
    )

	internal_nodes <- data.frame(
		name = add_node_names(phy)$node.label
	)

    if (!is.null(phy$node.label)) {
        internal_nodes$label = phy$node.label
    }
    
    nodes <- bind_rows(tip_nodes, internal_nodes)
    edges <- data.frame(phy$edge)
    colnames(edges) <- c("from", "to")
    
    gtree <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
    
    return(gtree)
}
