library(shiny)
library(igraph)


edge_df <- read.table(
    "https://bergplace.org/share/out.radoslaw_email_email",
    header = FALSE,
    skip   = 2
  )[, 1:2]
colnames(edge_df) <- c("from", "to")

# Wszystkie ID węzłów
nodes <- sort(unique(c(edge_df$from, edge_df$to)))

edge_df$from <- factor(edge_df$from, levels = nodes)
edge_df$to   <- factor(edge_df$to,   levels = nodes)

# MACIERZ liczby maili cnt_ij: wiersze = nadawca, kolumny = odbiorca
cnt_mat <- xtabs(~ from + to, data = edge_df)
cnt_mat <- as.matrix(cnt_mat)   # matrix, nie table

# Łączna liczba maili wysłanych przez każdy węzeł i (cnt_i)
cnt_i <- rowSums(cnt_mat)

# MACIERZ wag wij = cnt_ij / cnt_i (wierszowo)
prob_mat <- matrix(0, nrow = nrow(cnt_mat), ncol = ncol(cnt_mat),
                   dimnames = dimnames(cnt_mat))
nonzero_rows <- which(cnt_i > 0)
prob_mat[nonzero_rows, ] <- cnt_mat[nonzero_rows, , drop = FALSE] /
  cnt_i[nonzero_rows]

g <- graph_from_adjacency_matrix(
  prob_mat,
  mode     = "directed",
  weighted = TRUE,
  diag     = FALSE
)

g <- simplify(
  g,
  remove.multiple = TRUE,
  remove.loops    = TRUE,
  edge.attr.comb  = list(weight = "first")
)

stopifnot(vcount(g) == 167, ecount(g) == 5783)

n_nodes   <- vcount(g)
seed_frac <- 0.05
n_seeds   <- ceiling(seed_frac * n_nodes)

outdeg <- degree(g, mode = "out")
btw    <- betweenness(g, directed = TRUE)
clo    <- closeness(g, mode = "out")
pr     <- page.rank(g, directed = TRUE)$vector

seeds_outdeg <- order(outdeg, decreasing = TRUE)[1:n_seeds]

seeds_btw <- order(btw, decreasing = TRUE)[1:n_seeds]

seeds_clo <- order(clo, decreasing = TRUE)[1:n_seeds]

seeds_pr <- order(pr, decreasing = TRUE)[1:n_seeds]


simulate_ic_matrix <- function(prob_mat,
                               seeds,
                               max_iter = 10,
                               prob_multiplier = 1) {
  n <- nrow(prob_mat)
  active       <- rep(FALSE, n)
  active[seeds] <- TRUE
  
  newly_active <- active
  
  tried <- matrix(FALSE, nrow = n, ncol = n)
  
  result <- numeric(max_iter + 1)
  result[1] <- sum(active)
  
  for (t in 1:max_iter) {
    sources <- which(newly_active)
    if (length(sources) == 0) {
      if (t < max_iter + 1) {
        result[(t + 1):(max_iter + 1)] <- result[t]
      }
      break
    }
    
    newly_active <- rep(FALSE, n)
    
    for (u in sources) {
      p_row <- prob_mat[u, ] * prob_multiplier
      
      p_row[p_row > 1] <- 1
      
      # możemy próbować tylko:
      # - jeśli jest krawędź (p > 0)
      # - jeśli nie jest jeszcze aktywny
      # - jeśli nie próbowaliśmy i->j wcześniej
      candidates <- which(p_row > 0 & !active & !tried[u, ])
      
      if (length(candidates) == 0) next
      
      rand_vals <- runif(length(candidates))
      success   <- rand_vals < p_row[candidates]
      
      tried[u, candidates] <- TRUE
      
      if (any(success)) {
        new_nodes <- candidates[success]
        active[new_nodes]       <- TRUE
        newly_active[new_nodes] <- TRUE
      }
    }
    
    result[t + 1] <- sum(active)
  }
  
  result
}

run_diffusion_experiments <- function(prob_mat,
                                      max_iter    = 10,
                                      prob_slider = 100,
                                      n_runs      = 100,
                                      seed_frac   = 0.05) {
  prob_multiplier <- prob_slider / 100
  n <- nrow(prob_mat)
  n_seeds <- ceiling(seed_frac * n)
  
  iterations <- 0:max_iter
  n_iter <- length(iterations)
  
  strategies <- c("outdegree", "betweenness", "closeness", "random", "pagerank")
  results <- matrix(NA_real_, nrow = n_iter, ncol = length(strategies))
  colnames(results) <- strategies
  
  for (s in strategies) {
    sim_mat <- matrix(NA_real_, nrow = n_runs, ncol = n_iter)
    
    for (r in 1:n_runs) {
      if (s == "outdegree") {
        seeds <- seeds_outdeg
      } else if (s == "betweenness") {
        seeds <- seeds_btw
      } else if (s == "closeness") {
        seeds <- seeds_clo
      } else if (s == "pagerank") {
        seeds <- seeds_pr
      } else if (s == "random") {
        seeds <- sample(seq_len(n), n_seeds)
      } else {
        stop("Unknown strategy: ", s)
      }
      
      sim_vec <- simulate_ic_matrix(
        prob_mat        = prob_mat,
        seeds           = seeds,
        max_iter        = max_iter,
        prob_multiplier = prob_multiplier
      )
      
      sim_mat[r, ] <- sim_vec
    }
    
    results[, s] <- colMeans(sim_mat)
  }
  
  list(
    iterations = iterations,
    avg_active = results
  )
}


ui <- fluidPage(
  titlePanel("Dyfuzja informacji w sieci e-mail"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        inputId = "prob_mult",
        label   = "Mnożnik prawdopodobieństwa aktywacji",
        min     = 10,
        max     = 200,
        value   = 100,
        step    = 10
      ),
      sliderInput(
        inputId = "max_iter",
        label   = "Liczba iteracji procesu",
        min     = 1,
        max     = 50,
        value   = 10,
        step    = 1
      ),
    ),
    
    mainPanel(
      plotOutput("diffusionPlot"),
      verbatimTextOutput("graphSummary")
    )
  )
)


server <- function(input, output, session) {
  
  output$graphSummary <- renderPrint({
    cat("Liczba węzłów:", vcount(g), "\n")
    cat("Liczba krawędzi:", ecount(g), "\n\n")
    cat("Liczba węzłów startowych (5%):", n_seeds, "\n")
  })
  
  output$diffusionPlot <- renderPlot({
    max_iter    <- input$max_iter
    prob_slider <- input$prob_mult
    
    set.seed(123)  # powtarzalność wyników
    
    res <- run_diffusion_experiments(
      prob_mat    = prob_mat,
      max_iter    = max_iter,
      prob_slider = prob_slider,
      n_runs      = 100,
      seed_frac   = seed_frac
    )
    
    iters <- res$iterations
    avg   <- res$avg_active
    
    matplot(
      x    = iters,
      y    = avg,
      type = "l",
      lty  = 1,
      lwd  = 2,
      xlab = "Numer iteracji",
      ylab = "Średnia liczba aktywnych węzłów",
      main = "Dyfuzja informacji dla różnych zestawów węzłów początkowych"
    )
    
    legend(
      "bottomright",
      legend = c(
        "Najwyższy outdegree",
        "Najwyższe betweenness",
        "Najwyższe closeness",
        "Losowe węzły",
        "Najwyższy PageRank"
      ),
      col  = seq_len(ncol(avg)),
      lty  = 1,
      lwd  = 2,
      bty  = "n"
    )
  })
}


shinyApp(ui = ui, server = server)
