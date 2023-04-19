library("testthat")
test_that(
  "Testing recode_traits_gen()",
  {

    td <- readRDS("../testdata/td1.rds")
    amal.deps <- readRDS("../testdata/amaldeps1.rds")


    td2 <- recode_traits_gen(td, amal.deps, tax.col = TRUE, as.is = TRUE)


    expect_true(all(sort(unique(names(td$dat)[-1])) == sort(unique(names(td2$dat)))))

    expect_true(all(class(td2) == c("treedata", "list")))

  }
)
