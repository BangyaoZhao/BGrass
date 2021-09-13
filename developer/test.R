i = 2
load(paste0('../simulation/data/simulation6/data', i, '.Rdata'))
load(paste0('../simulation/data/simulation6/true_values.Rdata'))
devtools::load_all()
###########################
BV_chain = new_BV_chain(
  A,
  V,
  X,
  L,
  aggregation = TRUE,
  enrichment = FALSE,
  pi_delta = 0.5
)
BV_chain1 = update.BV_chain(BV_chain, 100,50)
BV_chain2 = update.BV_chain(BV_chain, 100,50)

plot(logLik(BV_chain1))

sub_BV_chain(BV_chain1, 75,101)

BV_chain1$start_pos
BV_chain1$end_pos

BV_chains=as.BV_chains(list(BV_chain1,BV_chain2))

BV_summ=summary(BV_chains)

coef_lst = coef(BV_chains)
coef_lst$beta_marginalized

BV_chains_summ=summary(BV_chains)


mean((coef_lst$beta_marginalized-beta*delta)^2)
mean((BV_summ$beta_table[,1]-beta*delta)^2)
###########################
BV_chain_enrichment = new_BV_chain(
  A,
  V,
  X,
  L,
  G = G,
  enrichment = TRUE
)
BV_chain_enrichment1 = update.BV_chain(BV_chain_enrichment, 1000)
BV_chain_enrichment2 = update.BV_chain(BV_chain_enrichment, 1000)

BV_chains=as.BV_chains(list(BV_chain_enrichment1,
                            BV_chain_enrichment2))

