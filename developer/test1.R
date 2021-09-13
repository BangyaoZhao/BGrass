data.frame(
  betaBV = coef_lst$beta,
  betaBV_enrichment = coef_lst_enrichment$beta,
  beta_true = beta * delta,
  deltaBV = coef_lst$delta,
  deltaBV_enrichment = coef_lst_enrichment$delta
) ->
  data

data %>%
  kableExtra::kbl(digits = 3)


