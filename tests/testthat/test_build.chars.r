library("testthat")
test_that(
  "Testing build.chars()",
  {
    
    extchar <- readRDS("../testdata/extchar.rds")
    char_info <- readRDS("../testdata/char_info.rds")
    
    bdchar <- build.chars(extchar, char_info)
    
    expect_true(length(extchar$phenotypes) == (length(bdchar$solved$chars) + length(bdchar$unsolved$chars)))
    expect_true(length(extchar$phenotypes) == (length(bdchar$solved$clusters) + length(bdchar$unsolved$clusters)))
    expect_true(length(extchar$phenotypes) == (length(bdchar$solved$tokens) + length(bdchar$unsolved$tokens)))
    
  }
)