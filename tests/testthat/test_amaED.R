library("testthat")
test_that(
  "Testing amaED()",
  {
    Tl <- matrix(c(-1, 1, 1, -1), 2, 2,
                 byrow = TRUE,
                 dimnames =list( c("T*", "T"), c("T*", "T")))

    C <- matrix(c(-2, 2, 2, -2), 2, 2,
               byrow = TRUE,
               dimnames =list( c("r", "b"), c("r", "b")))

    out <- amaED(Tl, C, type=c("ql"), phi=NULL)

    expected_out <- rbind(
      c(-2, 1, 1),
      c(1, -3, 2),
      c(1, 2, -3)
    )

    rownames(expected_out) <- colnames(expected_out) <- c("T*", "Tr", "Tb")
    expect_identical(out, expected_out)

  }
)





