library("testthat")
test_that(
  "Testing is_slumpable()",
  {

   Q <- initQ(c(1, 2), c(.3,.2))
   Q.h <- EHEtransform(Q)

   part_scheme1 <- list(c(1, 2), c(3,4,5))

   expect_true(is_slumpable(Q.h, part_scheme1))

   Q.h <- CHEtransform(Q)$Q
   part_scheme3 <- list(c(1:4), c(5:8))
   expect_true(is_slumpable(Q.h, part_scheme3))

  }
)



