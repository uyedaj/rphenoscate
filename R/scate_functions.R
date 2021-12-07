

#####################################
##### General utility functions #####
#####################################


#' Recodes a treedata object based on amalgamated characters.
#'
#' @export
recode_td <- function(td, traits, states, depstates=numeric(0)){
  tmp <- select(td, traits)
  hidden0 <- names(depstates)
  obs_char <- which(!(1:length(traits) %in% hidden0))
  for(i in 1:ncol(tmp$dat)){
      parent <- tmp$dat[[depstates[traits[i]]]]
      if(is.null(parent)){
        missingObs <- "*"
        tmp$dat[[i]] <- recode(as.character(tmp$dat[[i]]), "1"="1", "0"="0", "1 and 0"=missingObs, "0 and 1"=missingObs, .missing=missingObs)
      } else {
        tmp2 <- as.character(tmp$dat[[i]])
        tmp2[parent=="0"] <- "*"
        tmp$dat[[i]] <- recode(tmp2, "1"="1", "0"="0", "1 and 0"="*", "0 and 1"="*", .missing="*")
      }
      #recode0 <- ifelse(i %in% hidden0, "*", "0")
      
      #tmp$dat[[i]] <- recode(as.character(tmp$dat[[i]]), "1"="1", "0"=recode0, "1 and 0"=missingObs, "0 and 1"=missingObs, .missing=missingObs)
  }
  new.char <- tidyr::unite(tmp$dat, "new", sep="")
  new.char <- unname(sapply(new.char[[1]], function(x) paste(which(grepl(glob2rx(x), states))-1, collapse="&")))
  new.char[which(new.char=="")] <- "?"
  
  new.td <- select(td, -one_of(traits))
  new.td$dat[[paste(traits, collapse="+")]] <- new.char
  return(new.td)
}


#' Amalgamates traits using the matrices returned by amalgamate_deps.
#' 
#' @param td A treeplyr treedata object with characters to recode 
#' @param M A list produced by the function `amalgamate_deps`
#' 
#' @export
recode_traits <- function(td, M){
  # iterate through each new amalgamated trait and recode the matrix for that trait into td
  for(i in seq_along(M$new_traits)){
    trait_name <- M$new_traits[[i]]
    
    # if trait is not connected (ie has no dependencies)
    if(pracma::strcmp(M$traits[[i]], M$new_traits[[i]]) == TRUE){
      states <- setNames(c("0", "1"), 1:2)
      td <- recode_td(td, trait_name, states)
      td$dat[[trait_name]][td$dat[[trait_name]]=="?"] <- paste0(((1:length(states))-1), collapse="&")
    }
    else{ # if it has dependencies, recode the combined matrix for that subgraph into td
      # set gtraits 
      gtraits <- M$traits[[i]]
      # set states
      states <- colnames(M$M[[trait_name]])
      depstates <- M$depstates[[trait_name]]
      td <-recode_td(td, gtraits, states, depstates) 
      td$dat[[trait_name]][td$dat[[trait_name]]=="?"] <- paste0(((1:length(states))-1), collapse="&")
    }
  }
  
  return(td)
}


#' Generates a graph and model for traits given a dependency matrix.
#' 
#' @param dep_mat A dependency matrix
#' 
#' @details The dependency matrix can come from `rphenoscape::pa_dep_matrix` or be user generated
#' 
#' @return A list of: 
#' $M Transition matrices for the Structured-Hidden State Markov Models that can be run in RevBayes, 
#' corHMM or other tools
#' $graph The full graph of all dependencies among traits
#' $con.comp A list of the group membership etc. produced by `igraph::components`
#' $traits The old trait names that were amalgamated
#' $new.traits The new names of the amalgamated trait names
#' $depstates A list providing the dependencies of all characters
#' @export
amalgamate_deps <- function(dep_mat) {
  M <- list() # store all of our combined matrices
  traits <- list() 
  new_traits <- list() # store new trait names
  
  # remove redundant dependencies
  dep_mat <- remove_indirect(dep_mat)
  
  # generate our graph
  g <- igraph::graph_from_adjacency_matrix(t(as.matrix(dep_mat)))
  con.comp <- igraph::components(g, "weak")
  hidden <- list()
  # break up our graph into subgraphs, and combine matrices for the traits in those graphs
  for(i in 1:con.comp$no){
    g_sub <- igraph::induced_subgraph(g, which(con.comp$membership == i))
    .M <- combine_subgraph(g_sub)
    #bin_mats <- init_binary_matrix(g_sub) 
    #gtraits <- names(bin_mats) # get names of traits
    #newTraitName <- paste(gtraits, collapse="+") # name of new, amalgamated trait
    newTraitName <- .M$newTraitName
    traits[[i]] <- .M$traits # store names of traits before they were amalgamated
    new_traits[[i]] <- newTraitName # store name of new trait
    
    hidden[[newTraitName]] <- .M$dep_state
    M[[newTraitName]] <-  .M$M # store new matrix
  }

  return(list(M=M, graph=g, con.comp=con.comp, traits=traits, new_traits=new_traits, depstates=hidden))
}


## Not exported ##
# combine binary matrices for a subgraph on a dependency graph
# return the final combined matrix of all traits
combine_subgraph <- function(subgraph){
  # sort graph so that the "root" comes first
  sorted_g <- igraph::topo_sort(subgraph, mode = "out")
  bin_mats <- init_binary_matrix(subgraph) 
  
  # setting up:
  M <- list() 
  node_list <- names(sorted_g) 
  # trait_order <- list()  # keeps track of our traits in order of their dependencies
  ancestor_v <- node_list[[1]] # set to root
  M[[1]] <- bin_mats[[ancestor_v]]
  # trait_order[[1]] <- ancestor_v
  # now combine matrices based on dependencies
  # for each node in subgraph, combine trait matrix with its ancestor matrix
  Deps <- get_parent_states(subgraph)
  for (i in seq_along(node_list)){
    if (i != 1) { # if we are on the root node don't combine
      node_name <- node_list[[i]]
    
      if(i == 2) { # if we are currently on the first non-root vertex, we want to combine our current binary matrix with the binary matrix of the root vertex
        # set dependent state
        dep_state <- get_dep_state(subgraph, igraph::V(subgraph)[name == node_name], bin_mats[[ancestor_v]], node_list) #trait_order)
        
        M[[i]] <- comb2matrices(bin_mats[[ancestor_v]], bin_mats[[node_name]], dependent.state = dep_state)
      }
      else { # otherwise, combine with the resulting matrix from the combination we just did
        # set dependent state
        dep_state <- get_dep_state(subgraph, igraph::V(subgraph)[name == node_name], M[[i - 1]], node_list) #trait_order)
        M[[i]] <- comb2matrices(M[[i - 1]], bin_mats[[node_name]], dependent.state = dep_state)
      }
      # for debugging:
      # print(M[[i]]) 
      # print(dep_state)
    }
  }
  newTraitName <- paste(node_list, collapse="+") # name of new, amalgamated trait
  M[[newTraitName]] <- M[[i]]
  return(list(M=M[[newTraitName]], ancestor=ancestor_v, traits=node_list, newTraitName=newTraitName, dep_state=Deps))
}


## Not exported ##
get_parent_states <- function(subgraph){
  # sort graph so that the "root" comes first
  sorted_g <- igraph::topo_sort(subgraph, mode = "out")
  node_list <- names(sorted_g) 
  adj_v <- lapply(node_list,function(x) igraph::adjacent_vertices(x, graph=subgraph)[[1]])
  name_v <- lapply(1:length(adj_v), function(x) rep(x, length(adj_v[[x]])))
  parents <- setNames(unlist(name_v), names(unlist(adj_v)))
  return(parents)
}


## Not exported ##
# given a graph and a vertex v on that graph, the trait matrix for the ancestor node of v,
# and a list containing the traits in order of their combination into the matrix, 
# return the columns in the ancestor trait matrix which v depends on
get_dep_state <- function(subgraph, v, ancestor_mat, trait_order) {
  # get colnames of ancestor matrix
  mat_cols <- colnames(ancestor_mat)
  
  ## find which position the dependent trait(s) occupies in the colnames of our trait matrix
  adj_v <- igraph::adjacent_vertices(subgraph, v, mode="in") 
  ancestor_v <- names(adj_v[[v$name]]) # store the list of ancestor vertices of v
  
  # create a list containing positions in the ancestor_mat character where we have dependencies
  trait_pos <- list() 
  for(i in seq_along(ancestor_v)){
    pos <- which(trait_order == ancestor_v[[i]]) 
    trait_pos[[i]] <- pos
  }
  
  # for each position in the ancestor_mat character that v depends on,
  # store the corresponding column and add to our list
  depstate_list <- list() 
  for(i in seq_along(trait_pos)){
    depstate <- substr(mat_cols, trait_pos[[i]], trait_pos[[i]]) # extract values at that position for each column 
    depstate_list[[i]] <- c(which(depstate == 1)) # return the columns where the value for our dependent trait is 1
  }
  depstate_list <- unlist(depstate_list)
  return(unique(depstate_list))
}


#' Remove indirect links in a dependency graph.
#' 
#' Function for removing redundant dependency edges from our graph
#' 
#' @param dependency_matrix A dependency matrix
#'  
#' @export
remove_indirect <- function(dependency_matrix){
  z <- dependency_matrix
  diag(z) <- NA
  tmp <-  z & (!logical_mult(z,z))
  tmp[is.na(tmp)] <- 1
  res <- z*tmp
  return(res)
}


## Not exported ##
logical_mult = function(x,y)
{
  I = dim(x)[1]
  K = dim(x)[2]
  J = dim(y)[2]
  
  z = rep(FALSE,I*J)
  dim(z) = c(I,J)
  for(i in 1:I)
  {
    for(j in 1:J)
    {
      for(k in 1:K)
      {
        z[i,j] = z[i,j] || (x[i,k] && y[k,j])
      }
    }
  }
  z
}


#############################################
##### Structured Markov Model functions #####
#############################################




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


#' A wrapper for corHMM that fits Structured Hidden State Models.
#' 
#' @param td A treeplyr treedata object
#' @param amalgamations A list produced by `amalgamate_deps`
#' @param ... Additional arguments passed to `corHMM::rayDISC`
#' 
#' @details This function fits the models contained in `amalgamations$M`
#' to each of the characters in `td` using the function `corHMM::rayDISC`. 
#' 
#' @return A list of fits for each character in the dataset
#' @export
amalgamated_fits_corHMM <- function(td, amalgamations, ...){
  trees <- list()
  fits <- list()
  for(i in 1:ncol(td$dat)){
    .M <- amalgamations$M[[i]]
    diag(.M) <- NA
    .M[.M==0] <- NA
    colnames(.M) <- rownames(.M) <- 0:(ncol(.M)-1)
    dat <- data.frame(td[,i,tip.label=TRUE])
    #.testdat <- corHMM:::factorData(dat)
    .phy <- phytools::bind.tip(td$phy, "...delete...", edge.length = 0, where=length(td$phy$tip.label)+1)
    new.row <- data.frame("...delete...", paste(0:(ncol(.M)-1), collapse="&"))
    colnames(new.row) <- colnames(dat)
    .dat <- rbind(dat, new.row)
    fits[[i]] <- corHMM::rayDISC(.phy, .dat, rate.mat = .M, ...)
  }
  attributes(fits)$td <- td
  names(fits) <- colnames(td$dat)
  return(fits)
}


#' A wrapper for phytools that draws stochastic character maps given a corHMM fit.
#' 
#' @param fits A list produced by `amalgamated_fits_corHMM`
#' @param ... Additional arguments passed to `phytools::make.simmap`
#' 
#' @details This function takes the models produced from `amalgamated_fits_corHMM` and
#' uses them to generate stochastic character maps. 
#' 
#' @return A list of stochastic character maps for each character in the dataset
#' 
#' @export
amalgamated_simmaps_phytools <- function(fits, ...){
  td <- attributes(fits)$td
  trees <- list()
  for(i in 1:ncol(td$dat)){
    xx <- fits[[i]]$tip.states[1:length(td$phy$tip.label),]
    root <- fits[[i]]$tip.states[length(td$phy$tip.label)+1,]
    rownames(xx) <- td$phy$tip.label
    Q <- fits[[i]]$solution
    Q[is.na(Q)] <- 0
    diag(Q) <- -1*apply(Q, 1, sum)
    trees[[i]] <- phytools::make.simmap(td$phy, xx, Q=Q, pi=root, ...)
  }
  return(trees)
}

