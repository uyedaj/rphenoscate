library("testthat")
test_that(
  "Testing CHEtransform()",
  {

    Q <- initQ(c(1, 2),
               c(.3,.2)
    )

    tr <- CHEtransform(Q)

    tr$Q

    expect_true(all(sort(unique(c(round(tr$Q, digits = 1)))) == c(-0.3, -0.2, 0.0, 0.1)))

    ns <- c("1.1",  "1.2" , "1.3",  "1.4",  "2.1" , "2.2",  "2.3",  "2.4")

    expect_true(all(rownames(tr$Q) == ns))
    expect_true(all(colnames(tr$Q) == ns))

  }
)
