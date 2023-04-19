library("testthat")
test_that(
  "Testing comb2matrices()",
  {

    M1<-matrix(c(-1,1,  2,-2),2,2,byrow=TRUE)
    rownames(M1)<-colnames(M1)<-c("0","1")
    M2<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
    rownames(M2)<-colnames(M2)<-c("0","1")

    t1 <- rbind(c(-4, 3, 1, 0),
                c(4, -5, 0, 1),
                c(2, 0, -5, 3),
                c(0, 2, 4, -6))
    colnames(t1) <- rownames(t1) <- c("00", "01", "10", "11")

    t2 <- rbind(c(-1, 0, 1, 0),
                c(0, -1, 0, 1),
                c(2, 0, -5, 3),
                c(0, 2, 4, -6))
    colnames(t2) <- rownames(t2) <- colnames(t1)

    testobj1 <- comb2matrices(M1, M2, dependent.state = NULL)
    testobj2 <- comb2matrices(M1, M2, dependent.state = 2)

    expect_identical(testobj1, t1)
    expect_identical(testobj2, t2)
  }
)





