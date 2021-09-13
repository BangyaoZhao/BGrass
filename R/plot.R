# library(tidyverse)
# library(ggplot2)
# data.frame(BV_summ$lik$trace,
#            n_ite = BV_summ$chain_range) %>%
#   rename_with(function(x) {
#     paste0('chain', substr(x, 2, 10))
#   }, starts_with('X')) %>%
#   pivot_longer(starts_with('chain'),
#                names_to = 'chain',
#                values_to = 'logl') %>%
#   ggplot(aes(n_ite, logl, color = chain)) +
#   geom_line() +
#   geom_hline(yintercept = BV_summ$lik$each_chain[1], color = 'red') +
#   geom_hline(yintercept = BV_summ$lik$each_chain[2], color = 'red') +
#   geom_hline(yintercept = BV_summ$lik$all_chains) +
#   theme_bw()
