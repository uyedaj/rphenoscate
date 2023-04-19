library("testthat")
test_that(
  "Testing remove_indirect()",
  {

    m <- diag(5)
    diag(m) <- NA
    m[lower.tri(m)] <- 1

    m2 <- remove_indirect(m)

    expect_true(all(diag(m2[2:5, 1:4]) == 1))
    expect_true(all(m2[2:5, 1:4][lower.tri(m2[2:5, 1:4])] == 0))

  }
)
