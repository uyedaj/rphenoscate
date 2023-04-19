library("testthat")
test_that(
  "Testing filter_coverage()",
  {

    td <- readRDS("../testdata/fctd.rds")

    expect_identical(td, filter_coverage(td, traits = 0, taxa = 0))

  }
)

