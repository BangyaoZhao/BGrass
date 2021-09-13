#' Get the Laplacian matrix
#'
#' @param AEs The names of all the adverse events
#' @param data A data.frame with two columns
#' \describe{
#'   \item{AE}{The name of the AE}
#'   \item{group}{The name of the AE group}
#' }
#'
#' @return The Laplacian matrix induced by the AE groups
#'
#' @importFrom igraph graph.empty add.edges simplify laplacian_matrix
#' @import dplyr
#'
#' @export getLaplacian

getLaplacian = function(AEs, data) {
  # browser()
  J = length(AEs)

  data %>%
    distinct(group, AE) %>%
    filter(AE %in% AEs) %>%
    rowwise() %>%
    mutate(AE_ID = which(AE == AEs)) %>%
    ungroup() ->
    data

  AE_groups = split(data$AE_ID, data$group)

  g = graph.empty(J, directed = FALSE)

  v = numeric()
  for (AE_group in AE_groups) {
    if (length(AE_group) >= 2) {
      v = c(v, as.vector(combn(AE_group, 2)))
    }
  }

  g = add.edges(g, v)
  g = simplify(g)
  L = laplacian_matrix(g, normalized = TRUE)
  rownames(L) = AEs
  colnames(L) = AEs

  L = as.matrix(L)
  return(L)
}
