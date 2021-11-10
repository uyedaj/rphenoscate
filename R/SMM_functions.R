
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
#' @export
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



#' @title Initialize binary matrices given graph
#' @description Call matrices are populated with different parameters
#' @param graph igraph object
#' @return List of matrices
#' @examples
#' init_binary_matrix(g)
#' @export
init_binary_matrix<-function(graph){
  matrix.list=list()
  n.matrix=igraph::vcount(graph)
  vertices=igraph::V(graph)$name
  
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
#' @details
#' $binary.matrices # intial binary matrices assigned to each node of graph
#' $comb.matrices$matrix # combined matrix for each node
#' $comb.matrices$state.string # vector of states [1] "00" "01" "10" "11"
#' $comb.matrices$state.ident # specifies the order of ontology terms in each state [1] "UBERON:0007829" "UBERON:2000663"
#' $comb.matrices$state.observable # ids and names of "observable" states. In red-blue tail notation refers to blue and red states
#' $comb.matrices$state.hidden # ids and names of "hidden" states. In red-blue tail notation refers to "blue absent"
#' and "red absent"
#' $nodes.sorted # topologically sorted nodes
#' $vertex.hier # hierrachy of the nodes
#' @return List of matrices and their attributes
#' @examples
#' get_graph_matrix(g)
#' @export
get_graph_matrix<-function(graph){
  
  g=graph  
  
  # dependent chars object
  complex.char<-list()
  complex.char$binary.matrices<-init_binary_matrix(g)
  complex.char$comb.matrices<-list()
  
  # traverse graph
  .topo <- igraph::topo_sort(g, mode = c("out"))
  topo= names(.topo)
  complex.char$nodes.sorted=topo
  vertex.hier=ego(g, order=1, nodes = topo, mode = c("in"), mindist = 1)
  names(vertex.hier)<-topo
  complex.char$vertex.hier=vertex.hier
  
  
  for (i in seq_along(topo)){
    focal.v=complex.char$vertex.hier[[i]]
    
    if (length(focal.v)==0){
      complex.char$comb.matrices[[topo[i]]]$matrix=complex.char$binary.matrices[[topo[i]]]
      complex.char$comb.matrices[[topo[i]]]$state.string=row.names(complex.char$binary.matrices[[topo[i]]])
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
      r.name <- row.names(cmb)
      complex.char$comb.matrices[[topo[i]]]$state.string<-r.name #<-cmb %>% row.names()
      
      st.iden=unlist(lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$state.ident))
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
} 



#' @title Combining multiple matrices
#' @description Liely as independently evolving
#' @param list.matrices Lsit of matrices
#' @return Matrix
#' @export
combNmatrices<-function(list.matrices,  ...){
  comb.matrix<-list.matrices[[1]]
  
  for (i in 1:(length(list.matrices)-1)){
    comb.matrix=comb2matrices(comb.matrix, list.matrices[[i+1]], dependent.state=NULL, diag.as = 0)
  }
  
  return(comb.matrix)
}



#' @title Initialize binary matrices
#' @param char.state names for character states
#' @param rate.param names for the rate parameters
#' @param diag.as values to pas to the main diagonal elements
#' @return matrix
#' @export

init_char_matrix<-function(char.state, rate.param, diag.as=NA){
  n.state<-length(char.state)
  Q=matrix(ncol = n.state, nrow=n.state,byrow=TRUE)
  Q[xor(lower.tri(Q, diag = FALSE), upper.tri(Q, diag = FALSE))]<-rate.param
  Q<-t(Q)
  diag(Q)<-diag.as
  rownames(Q)<-colnames(Q)<-as.character(char.state)
  return(Q)
}

#' @title Make rate matrix for RevBayes
#' @param M rate matrix
#' @param prior prior to be used
#' @return matrix
#' @examples
#' cat(Mk_Rev(M))
#' @export


Mk_Rev<-function(M, prior="dnExp(20)"){
  
  dec<-paste("
NUM_STATES=", nrow(M), "\n",
             "for (i in 1:NUM_STATES) {
  for (j in 1:NUM_STATES) {
      rates[i][j] <-0.0
  }
}\n", sep="", collapse="")
  
 
  #i=1
  
  rt.dist<-paste("r", which(M>0), "~", prior, "\n", sep="")
  rt.moves<-paste("moves[++mvi] = mvScale(r", which(M>0), ", lambda=1, tune=true, weight=2)", "\n", sep="")
  
  rt<-c()
  for (i in 1:nrow(M))
    {
    which(M[i,]>0)->IJ
    
    rates<-M[i,IJ]
    
    #if (eq.diag==TRUE) rates<-M[i,IJ]/sum(M[i,IJ])
    
    # convert integer to decimal (RevBAyes sometimes retuns error if decimals are not used)
    #rates[rates%%1==0]<-paste(rates[rates%%1==0], ".0", sep="")
    
    
    rt<-c(rt,
          paste("rates[", i, "][", IJ, "]:=", sep ="" )
    )
    }
  
  rt<-paste(rt, "r", which(M>0), "\n", sep="")
  
  rt<-paste(rt, collapse="")
  rt.dist<-paste(rt.dist, collapse="")
  rt.moves<-paste(rt.moves, collapse="")
  
  rt<-paste(dec, rt.dist, rt.moves, rt, sep = "")
  
  rt.build <- '\nrate_matrix := fnFreeK(transition_rates=rates, rescaled=false, matrixExponentialMethod="eigen")\n'
  rt.root <- paste('root_freq <- simplex(', paste(rep(1, nrow(M)), collapse=" ,"), ')\n', sep="")
  
  rt <- paste(rt, rt.build, rt.root, sep="")

  return(rt)
}
#########

