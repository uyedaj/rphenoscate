library("testthat")
test_that(
  "Testing recode_traits()",
  {

    td <- readRDS("../testdata/td3.rds")
    amal.deps <- readRDS("../testdata/amaldeps3.rds")

    rctraits <- recode_traits(td, amal.deps, as.is = TRUE)

    expect_true(all(attributes(rctraits)$names == c("phy", "dat")))

  }
)
