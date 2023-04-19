library("testthat")
test_that(
  "Testing EHEtransform()",
  {

    Q <- initQ(c(1, 2), c(.3,.2))
    Qehe <- EHEtransform(Q)
    Qehe <- round(Qehe, digits = 1)

    EQ <- rbind(
      c(-0.3,  0.0,  0.1,  0.1,  0.1),
      c( 0.0 ,-0.3,  0.1,  0.1,  0.1),
      c( 0.1,  0.1, -0.2,  0.0,  0.0),
      c( 0.1  ,0.1 , 0.0 ,-0.2,  0.0),
      c( 0.1,  0.1,  0.0,  0.0, -0.2)
    )

    expect_true(all(Qehe == EQ))

  }
)
