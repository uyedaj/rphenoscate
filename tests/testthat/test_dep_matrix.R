library("testthat")
test_that(
  "Testing dep_matrix()",
  {

    ont <- ontologyIndex::get_ontology("../testdata/HAO.obo", extract_tags = "everything", propagate_relationships = c("BFO:0000050", "is_a"))

    tree <- readRDS("../testdata/td1.rds")
	
	tree$dat <- tree$dat[,1:6]
	
	colnames(tree$dat) <- c("taxon", "mandible", "acetabular groove", "paraocular carina", "tentorium", "anterior tentorial arm")

    dm <- dep_matrix(tree, ont, tax.col = TRUE)

    expect_true(all(colnames(dm) == rownames(dm)))
    expect_true(all(na.omit(unique(c(dm)))[1:2] %in% c(0, 1)))

  }
)
