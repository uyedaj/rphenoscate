library("testthat")
test_that(
  "Testing get_phi()",
  {

    A <- matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("A*", "A"), c("A*", "A")) )

    expect_true(all(get_phi(colnames(A)) == c(1,0)))

  }
)






