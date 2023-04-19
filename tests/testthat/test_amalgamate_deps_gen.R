library("testthat")
test_that(
  "Testing amalgamate_deps_gen()",
  {


    td <- readRDS("../testdata/td1.rds")
    dep.mat <- readRDS("../testdata/depmat1.rds")
    state.data <- readRDS("../testdata/statedata.rds")

    amal.deps <- amalgamate_deps_gen(td, dep.mat, mode = "check", state.data = state.data)


    expect_true(all(names(amal.deps$M) %in% amal.deps$traits))

    expect_true(all(attributes(amal.deps)$names == c("traits", "drop", "groups", "M", "states", "state.data")))


  }
)
