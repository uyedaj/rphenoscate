library("testthat")
test_that(
  "Testing amaSMM()",
  {

    C1 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
    C2 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("o", "g"), c("o", "g")) )
    C3 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("w", "m"), c("w", "m")) )
    arm <- amaSMM(C1, C2, C3)


    proper <- rbind(
      c(-3, 1 , 1, 0 , 1,  0 , 0,  0),
      c( 1,-3, 0, 1 , 0,  1 , 0,  0),
      c( 1, 0 ,-3, 1 , 0,  0 , 1,  0),
      c( 0, 1 , 1,-3, 0,  0 , 0,  1),
      c( 1, 0 , 0, 0 ,-3,  1 , 1,  0),
      c( 0, 1 , 0, 0 , 1, -3, 0,  1),
      c( 0, 0 , 1, 0 , 1,  0 ,-3,  1),
      c( 0, 0 , 0, 1 , 0,  1 , 1, -3)
    )

    expect_true(all(proper == arm))

  }
)
