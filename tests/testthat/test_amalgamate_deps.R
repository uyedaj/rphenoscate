library("testthat")
test_that(
  "Testing amalgamate_deps()",
  {

    names <- LETTERS[1:4]
    depmat <- diag(4)
    rownames(depmat) <- colnames(depmat) <- names
    depmat[1,2] <- depmat[3,4] <- 1

    adm <- amalgamate_deps(depmat)

    expect_condition(all(unlist(adm$new_traits)) == c("B+A", "D+C"))


    expect_true(adm$depstates$`B+A` == 1)
    expect_true(adm$depstates$`D+C` == 1)


  }
)


