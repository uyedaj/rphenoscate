#
# FUNCTIONS FOR WORKING WITH MATRICES
#

library("copula", lib.loc="~/.local/R/site-library") #stirling numbers
setwd("/home/tarasov/Dropbox/Rev-Bayes-test")
library(R.basic)
library("magrittr", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")


#' @title Combining two matrices
#' @description Combining two matrices. The parametric schem of matrice is defined by nattural
#' numbers; different numbers = different rate parameters
#' @param M1 matrix; if dependency true thenM1 controls M2
#' @param M2 matrix; if dependency true then: M2 depends on those states of M1 specified in dependent.state
#' @param dependent.state state(s) of M1 that switch on matrix M2 
#' @param name.sep separator for state names
#' @param diag.as hpopulate main diagonal with
#' @return Matrix
#' @examples
#' M1<-matrix(c(-1,1,  2,-2),2,2,byrow=TRUE)
#' rownames(M1)<-colnames(M1)<-c("0","1")
#' M2<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
#' rownames(M2)<-colnames(M2)<-c("0","1")
#' comb2matrices(M1, M2, dependent.state=NULL)
#' comb2matrices(M1, M2, dependent.state=2)

# if dependency true then: M2 depends on M1 states specified in dependent.state
comb2matrices<-function(M1,M2, dependent.state=NULL, name.sep="", diag.as=""){
  
  if (!is.null(dependent.state)){
  matrix.diag<-rep(0, ncol(M1))
  matrix.diag[dependent.state]<-1
  matrix.diag<-diag(matrix.diag)
  }
  
  if (is.null(dependent.state)){
    matrix.diag<-diag(nrow(M1))
    }
  
  M_kr=(M1%x%diag(nrow(M2)))+ (matrix.diag%x%M2)
  
  
  #getting colnames
  
  col=paste(colnames(kronecker(M1, diag(nrow(M2)), make.dimnames = T)),
            colnames(kronecker(diag(nrow(M1)), M2, make.dimnames = T)), sep="")
  col=gsub("::", name.sep, col, fixed=T)
  
  # merging two names
  rownames(M_kr)<-colnames(M_kr)<-col
  if (diag.as!="") diag(M_kr)<-diag.as
  
  return(M_kr)
}



###############################
#Declare functions
#####
indepen2matrices<-function(M1,M2, name.sep="", diag.as=""){
M_kr=(M1%x%diag(nrow(M2)))+ (diag(nrow(M1))%x%M2)

#getting colnames

col=paste(colnames(kronecker(M1, diag(nrow(M2)), make.dimnames = T)),
      colnames(kronecker(diag(nrow(M1)), M2, make.dimnames = T)), sep="")
col=gsub("::", name.sep, col, fixed=T)

# merging two names
rownames(M_kr)<-colnames(M_kr)<-col
if (diag.as!="") diag(M_kr)<-diag.as
return(M_kr)
}

#indepen2matrices(indepen2matrices(M1, M2, name.sep=":"), M1, name.sep=":")

indepen2matrices(M1, M2,  name.sep=":", diag.as=0)

#nrow(M1)

####################################
# populate all rows of coevolving matrix with different rates
########
rate_diff_joint<-function(M_kr){
N=which(M_kr!=0)%>%length
M_kr[which(M_kr!=0)]<-c(1:N)
M_kr=t(M_kr)
return(M_kr)
}

rate_diff_joint(M1)

##create matrices
m_dim=4
q=c(1,2,3,4,5,6,7,8)

M1<-matrix(c(-1,1,  2,-2),2,2,byrow=TRUE)
rownames(M1)<-colnames(M1)<-c("0","1")

M2<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
rownames(M2)<-colnames(M2)<-c("a","b")

M3<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
rownames(M2)<-colnames(M2)<-c("a","b")
# make them coevolving
M_kr=rate_diff_joint(indepen2matrices(M1, M2,  name.sep=":"))




#M_kr=indepen2matrices(M1, M2,  name.sep=":", diag.as=0)
#M_kr=indepen2matrices(M1, M2,  name.sep=":")
#M_kr=rate_diff_joint(M_kr)






##################
#
# check strong lumpability
#
###########
part_scheme=list(c(1, 2), c(3), c(4))
part_scheme=list(c(1), c(2), c(3,4))
part_scheme=list(c(1,2), c(3,4))
is_strg_lumpable(M_kr, list(c(1,2), c(3,4)))

M_kr=rows2rate_matrix(M.temp, row.vec, dt_rates[,7])
###
is_strg_lumpable<-function(M_kr, part_scheme){

Nper.part=lapply(part_scheme, length)%>%unlist
stat2M=matrix(0, nrow=length(Nper.part), ncol=ncol(M_kr))
for (i in 1:nrow(stat2M))
  stat2M[i, part_scheme[[i]]]<-1

M_rows=M_kr %*% t(stat2M)
tru.vals=c()

for (i in 1:length(part_scheme)){
  for (j in 1:ncol(M_rows)){
    tru.vals=c(tru.vals, length(unique(M_rows[part_scheme[[i]],j]))<2)
    }
}
# if this is false then matrix is not lumpable
return(all(tru.vals))
}

M3=comb2matrices(M1, M2, dependent.state=2)

M_kr=indepen2matrices(M1, M2,  name.sep=":", diag.as=0)
M_kr=indepen2matrices(M1, M2,  name.sep=":")
is_strg_lumpable(M3, part_scheme)

##########
# check weak lumpability: a special case involving division of rows bys sums of nrows


matrix_diff<-function(M_kr, part_scheme, t=1){

  Nper.part=lapply(part_scheme, length)%>%unlist
  stat2M=matrix(0, nrow=length(Nper.part), ncol=ncol(M_kr))
  for (i in 1:nrow(stat2M))
    stat2M[i, part_scheme[[i]]]<-1

 normalizer=diag(1/ Nper.part)
 diff=(normalizer %*% stat2M %*% expm::expm(M_kr*t) %*% t(stat2M))-
   (
     expm::expm(
       (normalizer %*% stat2M %*% M_kr %*% t(stat2M))*t
                 )
    )
  return(diff)
}

matrix_diff(M_kr, list(c(1, 2), c(3,4)), 10)
format(matrix_diff(M_kr, list(c(1, 2), c(3,4)), 20) , scientific=F)
max(abs(matrix_diff(M_kr, list(c(1, 2), c(3,4)))))
## max error between limpable matrix difference
lump_max_error<-function(...){
dif=matrix_diff(...)
max(abs(dif))
}
######
lump_max_error(M_kr, list(c(1, 2), c(3,4)))
#############################



