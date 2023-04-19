library("testthat")
library(tidyverse)

test_that(
  "Testing print_coverage()",
  {

    dat <- read.csv("../testdata/dat.csv")
    out <- print_coverage(dat)

    expect_lt(max(out$average, na.rm = TRUE), expected = 1.01)
    expect_gt(min(out$average, na.rm = TRUE), expected = -0.01)


  }
)

