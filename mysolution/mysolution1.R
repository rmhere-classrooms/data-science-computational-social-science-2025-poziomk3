# install.packages("igraph")
library(igraph)

set.seed(29112025)  # Ustawienie losowości dla powtarzalności wyników


n <- 100
p <- 0.05
g <- erdos.renyi.game(n = n, p.or.m = p)

summary(g)
is.weighted(g)   # graf nie jest ważony


V(g)               # lista wierzchołków

E(g)             # lista krawędzi


min_weight <- 0.01
max_weight <- 1

E(g)$weight <- runif(ecount(g), min = min_weight, max = max_weight)

summary(g)
is.weighted(g) # graf jest ważony


degree(g)
hist(degree(g),
     main = "Histogram stopni węzłów",
     xlab = "Stopień węzła",
     ylab = "Liczba węzłów")



clusters(g)


pr <- page.rank(g)$vector

plot(g, vertex.size=pr*500,
     vertex.label=NA, edge.arrow.size=.2)

