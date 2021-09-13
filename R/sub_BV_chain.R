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
  list2env(BV_chain, envir = environment())
  BV_chain$chain =
    resize_chain(BV_chain$chain,
                 new_start_pos - start_pos,
                 new_end_pos - end_pos)
  BV_chain$start_pos = new_start_pos
  BV_chain$end_pos = new_end_pos
  return(BV_chain)
}
