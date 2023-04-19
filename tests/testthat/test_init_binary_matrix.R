library("testthat")
test_that(
  "Testing init_binary_matrix()",
  {

    el <- cbind( c("jen", "tab", "joe"),
                 c("tab", "joe", "jen"))
    g <- graph.edgelist(el, directed=TRUE)
    bm <- init_binary_matrix(g)


    jen <- rbind(c(-1, 1),
                 c(2, -2))
    tab <- rbind(c(-3, 3),
                 c(4, -4))
    joe <- rbind(c(-5, 5),
                 c(6, -6))

    matches <- all((bm$jen == jen) &
                   (bm$tab == tab) &
                   (bm$joe == joe)
                   )

    expect_true(matches)
    expect_true(all(rownames(bm$jen) == c("0", "1")))
    expect_true(all(sapply(bm, rownames) == sapply(bm, colnames)))

  }
)

