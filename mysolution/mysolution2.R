# install.packages("igraph")  
library(igraph)

set.seed(29112025)  # dla powtarzalności


n <- 1000     

g <- barabasi.game(n = n,
                   directed = FALSE)

cat("=== Podsumowanie grafu BA ===\n")
print(summary(g))



plot(
  g,
  layout       = layout.fruchterman.reingold(g),
  vertex.size  = 3,
  vertex.label = NA,
)



btw <- betweenness(g)

max_btw <- max(btw)
which(btw == max_btw) # centralne węzeł (numer)
max_btw # wartość betweenness



diameter(g, directed = FALSE, weights = NA) # wartość średnicy grafu


###  Różnice między Barabási–Albert a Erdős–Rényi
# Graf Erdős–Rényi (ER):
# - Każda para węzłów jest łączona krawędzią z jednakowym prawdopodobieństwem p,
#   niezależnie od stopni węzłów.
#   Większość węzłów ma podobny, "średni" stopień.
# - Brak wyraźnych "hubów" – sieć jest bardziej jednorodna.

# Graf Barabási–Albert (BA):
# - Sieć rośnie w czasie: nowe węzły dołączają do już istniejących.
# - Nowy węzeł łączy się chętniej z węzłami,
#   które już mają wysoki stopień
# - Rozkład stopni jest zbliżony do potęgowego,
#   co oznacza obecność kilku węzłów-hubów o bardzo dużym stopniu
#   oraz wielu węzłów o małym stopniu.
# - BA lepiej modeluje wiele realnych sieci (np. sieci społeczne, WWW),
#   gdzie istnieją węzły o bardzo dużej liczbie połączeń.
