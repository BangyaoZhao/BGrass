library(pgdraw)
library(BayesLogit)
b_long <- readRDS("data/b_long.rds")
c_long <- readRDS("data/c_long.rds")
b <- readRDS("data/b.rds")
c <- readRDS("data/c.rds")

t1 = Sys.time()
for (i in 1:10) {
  Omega = pgdraw(b_long, c_long)
}
t2 = Sys.time()
t2 - t1

t1 = Sys.time()
for (i in 1:10) {
  Omega = rpg.devroye(length(b_long),b_long,c_long)
}
t2 = Sys.time()
t2 - t1



t1 = Sys.time()
OmegaMat1 = matrix(0, length(b), 10)
for (i in 1:10) {
  OmegaMat1[, i] = pgdraw(b, c)
}
t2 = Sys.time()
t2 - t1

t1 = Sys.time()
OmegaMat2 = matrix(0, length(b), 10)
for (i in 1:10) {
  OmegaMat2[, i] = rpg(length(b), b, c)
}
t2 = Sys.time()
t2 - t1

apply(OmegaMat1,2,mean)
apply(OmegaMat2,2,mean)
apply(OmegaMat1,2,sd)
apply(OmegaMat2,2,sd)


