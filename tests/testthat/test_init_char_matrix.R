library("testthat")
test_that(
  "Testing init_char_matrix()",
  {

    test <- init_char_matrix(c("0", "1"), c("0", "1"), diag.as = NA)

    expected <- rbind(c(NA, "0"),
                      c("1", NA))

    colnames(expected) <- rownames(expected) <- c("0", "1")

    expect_identical(test, expected)

  }
)
