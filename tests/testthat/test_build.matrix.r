library("testthat")
test_that(
  "Testing build.matrix()",
  {
    
    char_mat <- readRDS("../testdata/char_mat.rds")
    pheno_obj <- readRDS("../testdata/pheno_obj.rds")
    char_info <- readRDS("../testdata/char_info.rds")
    bdchar <- readRDS("../testdata/bdchar.rds")
    
    tax <- unique(row.names(char_mat))
    
    inf_mat <- build.matrix(tax, pheno_obj, char_info, bdchar)
    
    tokens <- apply(inf_mat, 2, unique)
    tokens <- lapply(tokens, function(x) unique(unlist(stringr::str_extract_all(x, "\\d"))) )
    
    expect_true(dim(inf_mat)[1] == length(tax))
    expect_true(dim(inf_mat)[2] == length(bdchar$solved$chars))
    expect_true(all(mapply(x = bdchar$solved$tokens, y = tokens, function(x,y) all(x %in% y) )))
    
  }
)