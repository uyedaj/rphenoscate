library("testthat")

test_that(
  "Testing process.ontofast()",
  {

    skip(message = "This is an interactive function that cannot be appropriately unit tested.")

    nex <- readRDS("../testdata/matrix.rds")

    mat <- readRDS("../testdata/mat.rds")

    HAO <- get_OBO("../testdata/HAO.obo", extract_tags = "everything",
                   propagate_relationships = c("BFO:0000050", "is_a"))

    u <- readRDS("../testdata/pof_u.rdf")

    s.terms <- c("margin", "edge","ridge", "carina", "spine", "spur", "angle",
                 "depression", "sclerite", "tergite", "sternite")
    g.terms <- c("length", "distance", "diameter")


    z <- process.ontofast(u, HAO, s.terms = s.terms, g.terms = g.terms)

    expect_true(all(names(u) == names(z)))

  }
)
