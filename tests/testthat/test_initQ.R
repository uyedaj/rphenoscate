library("testthat")
test_that(
  "Testing initQ()",
  {

    test <- initQ(c('a', 'p'), c(1,2), diag.as = NA)
    expected <- rbind(c(NA, 1),
                      c(2, NA))
    rownames(expected) <- colnames(expected) <- c("a", "p")

    expect_identical(test, expected)

  }
)
