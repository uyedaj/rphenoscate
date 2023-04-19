library("testthat")
test_that(
  "Testing as_matrixRB()",
  {

    test_val_1 <- "[[0.00, 2.00, 1.00, 0.00],\n [2.00, 0.00, 0.00, 1.00],\n [1.00, 0.00, 0.00, 2.00],\n [0.00, 1.00, 2.00, 0.00]]\n"
    test_val_2 <- "[[0.0, q[2], q[1], 0.0],\n [q[2], 0.0, 0.0, q[1]],\n [q[1], 0.0, 0.0, q[2]],\n [0.0, q[1], q[2], 0.0]]\n"

    Tl <- matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("T*", "T"), c("T*", "T")) )
    C <-matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
    Q <- amaSMM(Tl, C)

    rb1 <- as_matrixRB(Q)
    rb2 <- as_matrixRB(Q, symb='q')

    expect_equal(rb1, test_val_1)
    expect_equal(rb2, test_val_2)

  }
)

