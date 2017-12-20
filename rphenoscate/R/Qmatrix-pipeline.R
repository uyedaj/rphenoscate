library("magrittr")
library("igraph")

##############################################################################################
# FUNCTIONS
##############################################################################################

#' @title Combining two matrices
#' @description Combining two matrices. The parametric schem of matrice is defined by nattural
#' numbers; different numbers = different rate parameters
#' @param M1 matrix; if dependency true thenM1 controls M2
#' @param M2 matrix; if dependency true then: M2 depends on those states of M1 specified in dependent.state
#' @param dependent.state state(s) of M1 that switches on matrix M2 
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



#' @title Combining multiple matrices
#' @description Liely as independently evolving
#' @param list.matrices Lsit of matrices
#' @return Matrix
#' @examples
combNmatrices<-function(list.matrices,  ...){
  comb.matrix<-list.matrices[[1]]
  
  for (i in 1:(length(list.matrices)-1)){
    comb.matrix=comb2matrices(comb.matrix, list.matrices[[i+1]], dependent.state=NULL, diag.as = 0)
  }
  
  return(comb.matrix)
}

#' @title Initialize binary matrices given graph
#' @description Call matrices are populated with different parameters
#' @param graph igraph object
#' @return List of matrices
#' @examples
#' init_binary_matrix(g)
init_binary_matrix<-function(graph){
  matrix.list=list()
  n.matrix=vcount(graph)
  vertices=V(graph)$name
  
  param.scheme<-matrix(seq(1:(2* n.matrix )), ncol = 2, byrow=TRUE)
  for (i in 1:n.matrix){
    matrix.list[[vertices[i]]]<-matrix(c(-param.scheme[i,1],param.scheme[i,1], 
                                         param.scheme[i,2],-param.scheme[i,2]),2,2,byrow=TRUE)
    rownames(matrix.list[[vertices[i]]])<-colnames(matrix.list[[vertices[i]]])<-c("0","1")
  }
  
  return(matrix.list)
}


#' @title Get all dependency matrices given a dependecy graph
#' @description Construct dependency matrices and their correponding attributes
#' @param graph igraph object of ontology terms
#' @return List of matrices and their attributes
#' @examples
#' get_graph_matrix(g)

########## Structure of the List
$binary.matrices # intial binary matrices assigned to each node of graph

$comb.matrices$matrix # combined matrix for each node
$comb.matrices$state.string # vector of states [1] "00" "01" "10" "11"
$comb.matrices$state.ident # specifies the order of ontology terms in each state [1] "UBERON:0007829" "UBERON:2000663"
$comb.matrices$state.observable # ids and names of "observable" states. In red-blue tail notation refers to blue and red states
$comb.matrices$state.hidden # ids and names of "hidden" states. In red-blue tail notation refers to "blue absent"
# and "red absent"

$nodes.sorted # topologically sorted nodes
$vertex.hier # hierrachy of the nodes
###########

get_graph_matrix<-function(graph){
  
  g=graph  
  
  # dependent chars object
  complex.char<-list()
  complex.char$binary.matrices<-init_binary_matrix(g)
  complex.char$comb.matrices<-list()
  
  # traverse graph
  topo=topo_sort(g, mode = c("out")) %>% names()
  complex.char$nodes.sorted=topo
  vertex.hier=ego(g, order=1, nodes = topo, mode = c("in"), mindist = 1)
  names(vertex.hier)<-topo
  complex.char$vertex.hier=vertex.hier
  
  
  for (i in seq_along(topo)){
    focal.v=complex.char$vertex.hier[[i]]
    
    if (length(focal.v)==0){
      complex.char$comb.matrices[[topo[i]]]$matrix=complex.char$binary.matrices[[topo[i]]]
      complex.char$comb.matrices[[topo[i]]]$state.string=complex.char$binary.matrices[[topo[i]]] %>% row.names()
      complex.char$comb.matrices[[topo[i]]]$state.ident=topo[i]
      #complex.char$comb.matrices[[topo[i]]]$dependency.true=2
      complex.char$comb.matrices[[topo[i]]]$state.observable=integer(0)
      complex.char$comb.matrices[[topo[i]]]$state.hidden=integer(0)
    }
    
    if (length(focal.v)>0){
      
      if (length(focal.v)==1){ # if length =1 the dependency is chain like
        MC=complex.char$comb.matrices[[names(focal.v)]]$matrix
        M=complex.char$binary.matrices[[topo[i]]]
        dps=ncol(MC)
        cmb=comb2matrices(MC, M, dependent.state=dps, diag.as = 0)
      }
      
      
      if (length(focal.v)>1){
        # sequentially combine multiple matrices as independently coevolving
        list.matrices=lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$matrix)
        #names(list.matrices)=names(focal.v)
        #list.matrices=list.matrices[1:2]
        comb.mt=combNmatrices(list.matrices)
        
        # combine matrices  from above with the focal node matrix;
        # the state bearing dependency is where all entities=1, i.e. the last state
        M=complex.char$binary.matrices[[topo[i]]]
        dps=ncol(comb.mt)
        cmb=comb2matrices(comb.mt, M, dependent.state=dps, diag.as = 0)
      }
      
      # adding attributes
      complex.char$comb.matrices[[topo[i]]]$matrix=cmb
      complex.char$comb.matrices[[topo[i]]]$state.string<-r.name<-cmb %>% row.names()
      
      st.iden=lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$state.ident) %>% unlist()
      complex.char$comb.matrices[[topo[i]]]$state.ident<-st.iden<-c(st.iden, topo[i])
      
      # get observable and "hidden" states of the focal node
      ln=length(st.iden)-1 # observable a/p are only those which have all states from other chrs=1
      obs=which(substr(r.name, 1, ln)==paste(rep(1, ln), collapse=""))
      names(obs)=r.name[obs]
      complex.char$comb.matrices[[topo[i]]]$state.observable=obs
      
      hid=(1:length(r.name))[-obs]
      names(hid)<-r.name[hid]
      complex.char$comb.matrices[[topo[i]]]$state.hidden=hid
      
      
      
    } #end if (length(focal.v)>0)
  } #end all
  
  return(complex.char)
} # end function


# END FUNCTIONS
################################################################################################################################

#########################
# 
# TUTORIAL
#
#########################

# read dependency table
chrs.uber<-read.csv("dependencies.txt", header=T,  stringsAsFactors = F, na.strings = "")
# convert the table to matrix
chrs.uber<-as.matrix(chrs.uber[,c(1,2)])

# make igraph object
g=graph_from_edgelist(chrs.uber)
# plot all dependecies
plot(g, vertex.color="green", vertex.size=1,
     vertex.label.dist=0.5, vertex.label.cex=0.5, vertex.label=NULL, edge.arrow.size=.5 )

# now let's select some components to get matrices for them
# there are a few ways to do it:

# 1. working with connected components
con.comp=components(g, "weak")
# let's get a component with three nodes
com.id=which(con.comp$csize==3)
comp=con.comp$membership[con.comp$membership==com.id] %>% names()
# crating a new graph for this component
g1=subgraph(g, comp)
# plotting it
plot(g1)

# getting all matrices for g1
mt=get_graph_matrix(g1)
str(mt)
# matrix for a term
mt$comb.matrices$`UBERON:2002076`$matrix



# 2. working with subgraphs
# let's have a look on the subgraph of 16 nodes
com.id=which(con.comp$csize==16)
comp=con.comp$membership[con.comp$membership==com.id] %>% names()
# craeting a new graph for this component
g2=subgraph(g, comp)
plot(g2)
# now I make a subgraph for "UBERON:2000663" that includes only the nearest neighbors
g3=make_ego_graph(g, order=1, nodes = "UBERON:2000663", mode = c("all"), mindist = 0)
g3=g3[[1]]
plot(g3)

# getting all matrices for g3
mt=get_graph_matrix(g3)

# get a list of all combined matrices
comb.mt=lapply(mt$nodes.sorted, function(x) mt$comb.matrices[[x]]$matrix)
names(comb.mt)=mt$nodes.sorted
comb.mt


#####################
#
# JUNK CODE
#
#####################

# faked toy exmple
chrs.uber<-matrix(c(
  "UBERON:2201587", "UBERON:2102027",
  "UBERON:2201587", "UBERON:4300092",
  "UBERON:2202028", "UBERON:2201587",
  "UBERON:4200103", "UBERON:2201587",
  "UBERON:4200105", "UBERON:2201587"
),
ncol=2, byrow = T
)



#sbc=subcomponent(g, "UBERON:2001537", mode="in")

# #sorting connected components and saving them
# ub.ids=cbind(names(con.comp$membership), as.numeric(unname(con.comp$membership)))
# ub.ids=ub.ids[order(ub.ids[,2]),]
# colnames(ub.ids)<-c("id", "component_category")
# 
# write.csv(ub.ids, file="chrs_as_connected_components.csv")
##########

# # sorting graph to navigate
# node.sorted=topo_sort(g, mode = c("out"))
# node.sorted=names(node.sorted)
# chrs.uber=chrs.uber[factor(chrs.uber[,1], levels=node.sorted) %>% order(),]
# # now graph form orderd table
# #g=graph_from_edgelist(chrs.uber)

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
