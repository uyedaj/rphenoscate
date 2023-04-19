library("testthat")
test_that(
  "Testing makeTraitTree()",
  {

    tree <- readRDS("../testdata/td1.rds")
    ont <- ontologyIndex::get_ontology("../testdata/HAO.obo", extract_tags = "everything", propagate_relationships = c("BFO:0000050", "is_a"))

    tt <- makeTraitTree(tree, external = TRUE, ONT = ont)


    exp_chars <- names(tree$dat)
    true_chars <- tt$tip.label

    expect_true(all(exp_chars == true_chars))
  }
)

