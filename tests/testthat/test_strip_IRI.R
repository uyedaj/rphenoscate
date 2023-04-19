library("testthat")
test_that(
  "Testing strip_IRI()",
  {

    IRIs <- c("http://purl.obolibrary.org/obo/UBERON_0002396","http://purl.obolibrary.org/obo/UBERON_2001930",
    "http://purl.obolibrary.org/obo/UBERON_0011634","http://purl.obolibrary.org/obo/UBERON_2002024",
    "http://purl.obolibrary.org/obo/UBERON_2001969","http://purl.obolibrary.org/obo/UBERON_2001968")

    out <- strip_IRI(IRIs)
    expected_forms <- stringr::str_detect(out, "^UBERON:[:digit:]+$")
    expect_true(all(expected_forms))

  }
)






