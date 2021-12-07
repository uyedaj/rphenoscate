#' Automatically recode amalgamated traits in PARAMO 
#' Amalgamates traits using the matrices returned by amalgamate_deps
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


#' Generates a graph and model for traits given a dependency matrix
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


get_parent_states <- function(subgraph){
  # sort graph so that the "root" comes first
  sorted_g <- igraph::topo_sort(subgraph, mode = "out")
  node_list <- names(sorted_g) 
  adj_v <- lapply(node_list,function(x) igraph::adjacent_vertices(x, graph=subgraph)[[1]])
  name_v <- lapply(1:length(adj_v), function(x) rep(x, length(adj_v[[x]])))
  parents <- setNames(unlist(name_v), names(unlist(adj_v)))
  return(parents)
}

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


#' Remove indirect links in a dependency graph
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

#' A wrapper for corHMM that fits Structured Hidden State Models
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

#' A wrapper for phytools that draws stochastic character maps given a corHMM fit
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


