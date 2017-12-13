library("magrittr", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
setwd("~/Documents/Phenoscate/Phenoscate-R/rphenoscate/rphenoscate/R")


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

######################################





# read dependecy table
chrs.depen<-read.csv("dependencies.txt", header=T,  stringsAsFactors = F, na.strings = "")

# get a dependent pair of chars
char.id<-2
seq.of.chars<-which(chrs.depen[,1]==char.id)
seq.of.chars<-chrs.depen[seq.of.chars,3]%>% unique()
# first chrs in the vecto depends on the following chars
seq.of.chars<-c(char.id, seq.of.chars)

#intialize N matrices where N=number of binary chars involved in a given dependecy scheme,
# all matrices are populated with different parameters
# this can be better done using array()
matrix.list=list()
#
# each row is the id of rate parameters which will be plugged in rate matrices
param.scheme<-matrix(seq(1:(2*length(seq.of.chars))), ncol = 2, byrow=TRUE)
for (i in seq_along(seq.of.chars)){
  matrix.list[[i]]<-matrix(c(-param.scheme[i,1],param.scheme[i,1], 
                             param.scheme[i,2],-param.scheme[i,2]),2,2,byrow=TRUE)
  rownames(matrix.list[[i]])<-colnames(matrix.list[[i]])<-c("0","1")
}

# combine matrices together given their dependecies
# this matrix can be used for inference with hidden models
combined.matrix=comb2matrices(matrix.list[[2]], matrix.list[[1]], dependent.state=2, diag.as=0)

# check lumpability in respect to given partitioning scheme for states reflecting meaningfull
# absence, so we can use non-hidden models
# I will work on this function to more to get it working properly
#part_scheme=list(c(1, 2), c(3), c(4))
#is_strg_lumpable(combined.matrix, part_scheme)
# if matrix is lumpable then we construct the aggregated chain, let's imagine the matrix
# is lumpable
is.lumpable=T

if (is.lumpable==T) {
  #identify which states are not meaningfull given absence
  # I will work on it to extend for arbitrary number of combined matrices
  # by far I just create a matrix how it should look like
  aggregated.matrix<-matrix(c(0,1,2,  3,0,4, 5, 6, 0),3,3,byrow=TRUE)
  rownames(aggregated.matrix)<-colnames(aggregated.matrix)<-c("0","10", "11")

}

### Outputs
combined.matrix # combined matrix
aggregated.matri # aggregated matrix
seq.of.chars # pair of chars where #1 depends on #2 (swith off dependency), also useful to 
#know character order in state names


#########################
# Topological sorting
#########################

g <- barabasi.game(100)
topo_sort(g)

el <- matrix( c("foo", "bar", "bar", "foobar"), nc = 2, byrow = TRUE)
g=graph_from_edgelist(el)
plot(g)

# make edge table without repeats
chrs.depen.red=chrs.depen[chrs.depen[,2]==1,]
chrs.edges<-as.matrix(chrs.depen.red[,c(3,1)])
rownames(chrs.edges)<-NULL
colnames(chrs.edges)<-NULL
as.character(chrs.edges)
chrs.edges<-matrix(as.character(chrs.edges), ncol = 2, byrow = F)

g=graph_from_edgelist(chrs.edges)
plot(g, vertex.color="green", vertex.size=1,
     vertex.label.dist=0.5, edge.arrow.size=.5 )

components(g, "weak")
######## datasets with uberon ids
chrs.uber<-read.csv("dependencies.txt", header=T,  stringsAsFactors = F, na.strings = "")
chrs.uber<-as.matrix(chrs.uber[,c(1,2)])
g=graph_from_edgelist(chrs.uber)
plot(g, vertex.color="green", vertex.size=1,
     vertex.label.dist=0.5, vertex.label=NA, edge.arrow.size=.5 )

con.comp=components(g, "weak")

#sorting connected components and saving them
ub.ids=cbind(names(con.comp$membership), as.numeric(unname(con.comp$membership)))
ub.ids=ub.ids[order(ub.ids[,2]),]
colnames(ub.ids)<-c("id", "component_category")

write.csv(ub.ids, file="chrs_as_connected_components.csv")
##########

# getting subgraph to navigate
topo_sort(g)

# get a subgraph for test
test=cbind(c(10, 10, 6)%>% as.character(), c(59, 60, 10)%>% as.character())
g=graph_from_edgelist(test)
plot(g)
node.sorted=topo_sort(g, mode = c("out"))
node.sorted=names(node.sorted)%>% as.numeric()
node.sorted=node.sorted[node.sorted%in%test[,1]]
match(node.sorted, test[,1])
test=cbind(test, match(test[,1], node.sorted))
test=test[order(test[,3]),]

# doing loop to merge matrices
g1=graph_from_edgelist(test[,c(1,2)])
plot(g1)

#make loop to get node ancestors
test1=test
i=1
#subcomp=matrix(NA, ncol=2, nrow=vcount(g1))
subcomp=list()
for (i in 1:vcount(g1)){
  tmp=subcomponent(g1, i, mode="in")%>% names()%>% as.numeric()
  subcomp[[i]]<-list(tmp[1], tmp[-1])
}

