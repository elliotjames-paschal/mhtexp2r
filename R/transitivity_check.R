#' Check for transitivity violations in hypothesis set
#'
#' @param subset_hypotheses List of new hypotheses
#' @param rejected List of previously rejected hypotheses
#'
#' @return TRUE if contradiction found, FALSE otherwise
#' @export
violates_transitivity <- function(subset_hypotheses, rejected) {
  all_hypotheses <- c(subset_hypotheses, rejected)

  # Group by (outcome, subgroup)
  keys <- vapply(all_hypotheses, function(h) {
    paste(h$outcome, h$subgroup, sep = "_")
  }, character(1))

  groups <- split(all_hypotheses, keys)

  for (group in groups) {
    # Convert to edge list, shift all treatment IDs by +1 to avoid 0s (igraph limitation)
    edges <- do.call(rbind, lapply(group, function(h) c(h$t1 + 1, h$t2 + 1)))

    # Construct undirected graph
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    clusters <- igraph::components(g)$membership

    # Check for logical contradiction: duplicate treatments in different clusters
    cluster_sets <- split(as.numeric(names(clusters)), clusters)
    treatments_flat <- unlist(cluster_sets)

    # If any treatment appears more than once, it’s in multiple clusters → contradiction
    if (length(unique(treatments_flat)) < length(treatments_flat)) {
      return(TRUE)
    }
  }

  return(FALSE)
}
