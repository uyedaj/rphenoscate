library("testthat")
test_that(
  "Testing Q2model()",
  {
    Q <- initQ(c(1, 2), c(.3,.2))
    test <- Q2model(Q)

    expected <- rbind(c(NA, 2),
                    c(1, NA))
    rownames(expected) <- colnames(expected) <- c("1", "2")

    expect_identical(test, expected)
  }
)

