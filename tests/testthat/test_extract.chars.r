library("testthat")
test_that(
  "Testing extract.chars()",
  {
    
    exclu <- readRDS("../testdata/exclu.rds")
    
    ch <- extract.chars(exclu)
    
    expect_true(dim(exclu$matrix)[1] == length(unlist(ch$phenotypes)))
    expect_true(dim(exclu$matrix)[1] == sum(unlist(lapply(ch$submatrices, function(x) sqrt(length(x)) ))))
    expect_true(dim(exclu$dataframe)[1] == ((length(unlist(ch$phenotypes))^2) - length(unlist(ch$phenotypes)))/2)
    
  }
)