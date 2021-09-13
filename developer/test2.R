devtools::load_all()

load("../Server/sim6/data/data1.Rdata")
load("../Server/sim6/data/true_values.Rdata")

BV_chain = new_BV_chain(A,
                        V,
                        X,
                        L,
                        enrichment = FALSE,
                        glm_start = TRUE)

BV_chain = update(BV_chain, n_ite = 20000)

l=logLik(BV_chain)
plot(l)


BV_coef=coef(BV_chain)
mean((BV_coef$beta-beta*delta)^2)


