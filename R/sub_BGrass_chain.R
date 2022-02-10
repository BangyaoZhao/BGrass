#' Extract a part of the BGrass chain
#'
#' @param BGrass_chain An object of s3 class BGrass_chain
#' @param new_start_pos New start position
#' @param new_end_pos New end position
#'
#' @return The truncated BGrass_chain
#'
#' @export sub_BGrass_chain

sub_BGrass_chain = function(BGrass_chain,
                        new_start_pos,
                        new_end_pos) {
  BGrass_chain$chain =
    BGrass_chain$chain[new_start_pos - start_pos + 1,
                   new_end_pos - start_pos + 1]
  BGrass_chain$start_pos = new_start_pos
  BGrass_chain$end_pos = new_end_pos
  return(BGrass_chain)
}
