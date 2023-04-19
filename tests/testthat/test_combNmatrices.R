library("testthat")
test_that(
  "Testing combNmatrices()",
  {

    M1<-matrix(c(-1,1,  2,-2),2,2,byrow=TRUE)
      rownames(M1)<-colnames(M1)<-c("0","1")
    M2<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
      rownames(M2)<-colnames(M2)<-c("0","1")
    M3<-matrix(c(-5,5,  6,-6),2,2,byrow=TRUE)
      rownames(M3)<-colnames(M3)<-c("0","1")

    ms <- list(M1=M1,
               M2=M2,
               M3=M3
          )

    combed <- combNmatrices(ms)


    expc <- rbind(
            c(0, 5, 3, 0, 1, 0, 0, 0),
            c(6, 0, 0, 3, 0, 1, 0, 0),
            c(4, 0, 0, 5, 0, 0, 1, 0),
            c(0, 4, 6, 0, 0, 0, 0, 1),
            c(2, 0, 0, 0, 0, 5, 3, 0),
            c(0, 2, 0, 0, 6, 0, 0, 3),
            c(0, 0, 2, 0, 4, 0, 0, 5),
            c(0, 0, 0, 2, 0, 4, 6, 0)
        )

    expect_true(all(combed == expc))

  }
)




