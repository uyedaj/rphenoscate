library("testthat")
test_that(
  "Testing auto_dep_matrix()",
  {

    tree <- ape::rtree(n = 12, tip.label = LETTERS[1:12])
    dats <- replicate(5,
                      rnorm(n=12, mean=sample(c(20, 40, 50, 60, 80)), sd = 15)
    )
    dats <- as.data.frame(dats)
    names(dats) <- letters[1:5]

    dats$Species <- tree$tip.label

    td <- suppressWarnings(treeplyr::make.treedata(tree, dats))

    auto_dep_matrix(td)

    expect_true(is.na(unique(diag(auto_dep_matrix(td)))))

    expect_true(unique(auto_dep_matrix(td)[upper.tri(auto_dep_matrix(td)) | lower.tri(auto_dep_matrix(td))]) == 0)

  }
)











quote({
  hi := 2
})













