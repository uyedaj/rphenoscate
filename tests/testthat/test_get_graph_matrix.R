library("testthat")
test_that(
  "Testing get_graph_matrix()",
  {

    g <- igraph::make_graph(c("A", "B", "B", "C", "C", "D"), directed = TRUE)
    m <- get_graph_matrix(g)


    check_format <- function(x){
      o <- all(diff(c(abs(x))) == c(1, -1, 1))
      t <- all(diag(x) < 0)
      return(o & t)
    }

    all_good <- all(sapply(m$binary.matrices, check_format))

    expect_true(all_good)

  }
)
