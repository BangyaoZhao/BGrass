#' Get the group matrix
#'
#' @param AEs The names of AEs
#' @param data A data frame that contains two columns
#' \itemize{
#'  \item{"AE"}{The names of the AE}
#'  \item{"group"}{The group that the AE belongs to}
#' }
#' @param min_size The minimal size of the group to be included in the matrix
#'
#' @export getG
#' @return The group matrix G

getG = function(AEs, data, min_size = 15) {
  # browser()
  J = length(AEs)
  groups = sort(unique(data$group))
  K = length(groups)
  G = matrix(0, J, K)
  rownames(G) = AEs
  colnames(G) = groups
  for (i in 1:nrow(data)) {
    G[data$AE[i], data$group[i]] = 1
  }
  G = G[, colSums(G) >= min_size]
  G = cbind(intercept = 1, G)
  return(G)
}
