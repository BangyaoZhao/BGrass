library(ggplot2)
library(dplyr)
i = 2
load(paste0('../simulation/data/simulation6/data', i, '.Rdata'))
load(paste0('../simulation/data/simulation6/true_values.Rdata'))
devtools::load_all()
###########################
BV_chain1 = new_BV_chain(
  A,
  V,
  X,
  L,
  aggregation = TRUE,
  enrichment = FALSE,
  pi_delta = 0.5,
  glm_start = TRUE
)
BV_chain1 = update.BV_chain(BV_chain1, 100)
###########################
BV_chain2 = new_BV_chain(
  A,
  V,
  X,
  L,
  aggregation = TRUE,
  enrichment = FALSE,
  pi_delta = 0.5,
  glm_start = FALSE
)
BV_chain2 = update.BV_chain(BV_chain2, 100)
###########################
l1 = logLik(BV_chain1, pos = 1:BV_chain1$chain_length)
l2 = logLik(BV_chain2, pos = 1:BV_chain2$chain_length)
data.frame(
  l = c(l1, l2),
  n_ite = rep(1:length(l1), 2),
  chain = paste0('chain', rep(1:2, each=length(l1)))
)%>%
  filter(n_ite>=0)%>%
  ggplot(aes(n_ite,l,color=chain))+
  #ylim(-1.6e5,NA)+
  geom_line()

chain_range=2000:10000
coef_lst1 = coef(BV_chain1,chain_range =chain_range )
coef_lst2 = coef(BV_chain2,chain_range = chain_range)
mean((coef_lst1$beta-beta*delta)^2)
mean((coef_lst2$beta-beta*delta)^2)

coef_lst2$tau_2

cbind(coef_lst1$beta,coef_lst2$beta)
