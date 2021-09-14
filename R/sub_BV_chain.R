#' Extract a part of the BV chain
#'
#' @param BV_chain An object of s3 class BV_chain
#' @param new_start_pos New start position
#' @param new_end_pos New end position
#'
#' @return The truncated BV_chain
#'
#' @export sub_BV_chain

sub_BV_chain = function(BV_chain,
                        new_start_pos,
                        new_end_pos) {
  BV_chain$chain =
    BV_chain$chain[new_start_pos - start_pos + 1,
                   new_end_pos - start_pos + 1]
  BV_chain$start_pos = new_start_pos
  BV_chain$end_pos = new_end_pos
  return(BV_chain)
}
