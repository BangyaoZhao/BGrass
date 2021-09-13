devtools::load_all()
chain_length = 5
chain = list(
  a_chain = 1:chain_length,
  beta_chain = matrix(rnorm(2 * chain_length), 2, chain_length),
  A_chain = array(rnorm(6 * chain_length), dim = c(3, 2, chain_length))
)

resize_chain(chain,1,0)$a_chain
