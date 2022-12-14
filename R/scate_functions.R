

#####################################
##### UPDATE: Diego (2022)      #####
#####################################

#' Obtains a presence-absence dependency matrix when using external ontologies.
#'
#' @param td A treeplyr treedata object with characters to recode
#' @param ONT An 'ontology_index' object of an external ontology imported in R using the package 'ontologyIndex'
#' @param tax.col Indicates if the first column of the data set contains taxon names
#'
#' @examples
#' \dontrun{
#'     dep_matrix(treedata, HAO, tax.col=TRUE)
#' }
#'
#' @import ontologyIndex
#'
#' @export
dep_matrix <- function(td, ONT, tax.col = FALSE){

  m <- as.data.frame(td$dat)

  if(tax.col == TRUE){ terms <- colnames(m)[-1] }else{ terms <- colnames(m) }

  iris <- names(ONT$name[match(terms, ONT$name)])

  z <- sapply(iris, function(x) intersection_with_descendants(ONT, roots = x, terms = iris) )

  mat <- do.call(cbind, lapply(z, function(x) iris %in% x ))
  mat <- mat*1
  diag(mat) <- NA

  colnames(mat) <- rownames(mat) <- terms

  return(mat)

}


#' Obtains a 'auto-dependency' matrix when amalgamating characters referring to the same anatomical entity.
#'
#' @param td A treeplyr treedata object with characters to recode
#' @param tax.col Indicates if the first column of the data set contains taxon names
#'
#' @examples
#' \dontrun{
#'     auto_dep_matrix(treedata, tax.col = FALSE)
#'  }
#'
#' @importFrom stats dist
#'
#' @export
auto_dep_matrix <- function(td, tax.col = FALSE){

  m <- as.data.frame(td$dat)

  if(tax.col == TRUE){
    mat <- as.matrix(dist(match(colnames(m)[-1], colnames(m)[-1]), diag = TRUE, upper = TRUE))
  }else{
    mat <- as.matrix(dist(match(colnames(m), colnames(m)), diag = TRUE, upper = TRUE))
    }

  mat[mat > 0] <- 2
  mat[mat == 0] <- 1
  mat[mat == 2] <- 0
  mat[upper.tri(mat)] <- 0
  diag(mat) <- NA

  if(tax.col == TRUE){
    colnames(mat) <- rownames(mat) <- colnames(m)[-1]
  }else{
    colnames(mat) <- rownames(mat) <- colnames(m)
    }

  return(mat)

}


#' Amalgamates and recodes traits using the matrices returned by amalgamate_deps_gen.
#'
#' @param td A treeplyr treedata object with characters to recode
#' @param amal.deps A list produced by the function `amalgamate_deps_gen`
#' @param tax.col Indicates if the first column of the data set contains taxon names
#'
#' @examples
#' \dontrun{
#'     recode_traits_gen(treedata, amal.deps)
#' }
#'
#' @export
recode_traits_gen <- function(td, amal.deps, tax.col = TRUE)
{

  # Format as data.frame #
  td$dat <- as.data.frame(td$dat)

  ## Prepare data ##
  # Get vectors for new traits #
  if(tax.col == T){ new.traits <- lapply(amal.deps$groups, function(x) { td$dat[, (x + 1)] } ) }else{

    new.traits <- lapply(amal.deps$groups, function(x) { td$dat[,x] } )

  }

  # Get index of non-dependent and dependent traits #
  nondep.index <- !sapply(new.traits, is.data.frame)
  dep.index <- sapply(new.traits, is.data.frame)

  ## Process non-dependent traits ##
  # Get non-dependent traits #
  nondep.traits <- new.traits[nondep.index]

  # Recode inapplicables for non-dependent traits #
  nondep.traits <- lapply(nondep.traits, function(x) gsub(x, pattern = "-", replacement = "?") )
  # OBS.: "-" (inapplicables) in non-dependent traits will be interpreted as "?" (missings) #

  # Get missing for non-dependent traits #
  miss.states <- lapply(amal.deps$states[nondep.index], function(x) paste0(x, collapse = "&") )

  # Recode missing for non-dependent traits #
  nondep.traits <- mapply(x = nondep.traits, y = miss.states, function(x,y){ gsub(x, pattern = "\\?", replacement = y) }, SIMPLIFY = F)

  ## Process dependent traits ##
  # OBS.: In this context 'dependency' means simply that two or more traits are annotated with the same ontology term.
  # For this more general type of dependency SMM models are used and traits recoded accordingly.
  # However, 'true' dependencies can also occur (i.e., parthood-relations or property-instantiation).
  # They require more complex models. In these cases, EDql and EDap models are used.

  # Get dependent traits #
  dep.traits <- new.traits[dep.index]

  # Merge states of dependent traits #
  dep.traits <- lapply(dep.traits, function(x) apply(x, 1, paste0, collapse = "") )

  # Get observed state combinations for dependent traits #
  obs.states <- lapply(dep.traits, function(x) unique(x) )

  # Get all possible state combinations for dependent traits #
  dep.states <- amal.deps$states[dep.index]

  # Get recoded states for dependent traits #
  recode.states <- lapply(dep.states, function(x) paste(0:(length(x) - 1)) )

  # Recoding of inapplicable states for dependent traits #
  for(i in 1:length(dep.traits)){

    if(length(grep(dep.traits[[i]], pattern = "0-")) > 0){

      dep.traits[[i]] <- gsub(dep.traits[[i]], pattern = "0-", replacement = "0")

    }

    if(length(grep(dep.traits[[i]], pattern = "1-")) > 0){

      dep.traits[[i]] <- gsub(dep.traits[[i]], pattern = "0\\1(.)", replacement = "1\\1")
      dep.traits[[i]] <- gsub(dep.traits[[i]], pattern = "1-", replacement = "0")

    }

  }

  # Recode all non-missings from dependent traits #
  new.states <- mapply(x = recode.states, y = dep.traits, z = dep.states, function(x,y,z) {x[match(y,z)]}, SIMPLIFY = F )

  # Get all elements with missing from depend traits #
  miss.states <- mapply(x = dep.traits, y = new.states, function(x,y) { x[is.na(y)] } )

  # Reorganize missing labels "?" #
  miss.states <- lapply(miss.states, function(x) { gsub(x, pattern = "\\?", replacement = ".") })

  # Get all recoded states from elements with missing from dependent traits #
  miss.states <- mapply(x = miss.states, y = recode.states, z = dep.states, function(x,y,z) { sapply(x, function(w) paste0(y[grepl(z, pattern = w)], collapse = "&") ) } )

  # Extract traits with missings from dependent traits #
  new.states.miss <- new.states[sapply(miss.states, is.character)]

  # Extract recoded states for traits with missing from dependent traits #
  miss.states <- miss.states[sapply(miss.states, is.character)]

  # Recode traits with missing from dependent traits #
  for(i in 1:length(new.states.miss)) { new.states.miss[[i]][is.na(new.states.miss[[i]])] <- miss.states[[i]] }

  # Join all recoded dependent traits #
  new.states[sapply(new.states, function(x) { any(is.na(x)) } )] <- new.states.miss

  # Join all recoded traits #
  new.traits[dep.index] <- new.states
  new.traits[nondep.index] <- nondep.traits

  # Build final data set #
  M <- as.data.frame(do.call(cbind, new.traits))
  td <- treeplyr::make.treedata(td$phy, cbind(td$dat[,1],M))

  # Return results #
  return(td)

}


#' Generates a model for traits given a dependency matrix when using external ontologies.
#'
#' @param td A treeplyr treedata object with characters to amalgamate
#' @param dep.mat A dependency matrix produced by the function `dep_matrix`
#' @param mode Indicates whether to perform automatic amalgamations using SMM models ('auto') or
#' to check for different types of dependencies using ED models ('check'). See Tarasov (2021) for details
#' @param state.data A list with information about character states descriptions. If mode = 'check', the function
#' will try to automatically check the type of trait dependency based of state descriptions
#' @param tax.col Indicates if the first column of the data set contains taxon names
#'
#' @examples
#'
#' \dontrun{
#'     amalgamate_deps_gen(treedata, dep_mat, state.data = state.data)
#' }
#'
#' @return A list of:
#' $traits The old trait names that were amalgamated
#' $drop The dependent traits that were dropped after amalgamation
#' $groups The groups of traits that were amalgamated
#' $M Transition matrices for the Markov Models (SMM or ED) that can be run in RevBayes, corHMM or other tools
#' $states A list providing the final amalgamated states
#' $state.data A list providing the character state descriptions for the amalgamated traits
#'
#' @export
amalgamate_deps_gen <- function(td, dep.mat, mode = c("auto", "check"), state.data, tax.col = TRUE)
{

  # Format as data.frame #
  td$dat <- as.data.frame(td$dat)

  # Create list to store results #
  amal.deps <- list()

  # Get all trait terms #
  if(tax.col == T){ amal.deps$traits <- colnames(td$dat)[-1] }else{ amal.deps$traits <- colnames(td$dat) }

  # Get total number of traits #
  N <- length(amal.deps$traits)

  # Get traits to drop based on dependencies #
  amal.deps$drop <- rowSums(dep.mat, na.rm = T)

  # Get groups of traits to keep #
  amal.deps$groups <- apply(dep.mat, 2, function(x) which(x != 0 | is.na(x)) )

  # Create several lists to store information #
  amal.deps$M <- list()
  amal.deps$states <- list()
  amal.deps$state.data <- list()

  ## Automatic DDA algorithm  across matrix (tentatively!) ##
  for(i in 1:N){

    # Get a single trait (or groups of traits to be amalgamated) #
    if(tax.col == T){ states <- td$dat[,(amal.deps$groups[[i]] + 1)] }else{ states <- td$dat[,amal.deps$groups[[i]]] }

    # Get state information #
    amal.deps$state.data[[i]] <- state.data[amal.deps$groups[[i]]]

    # Process groups of dependent traits #
    if(length(dim(states)) > 0){

      states <- lapply(apply(states, 2, unique, simplify = F), function(x) x[order(x)] )
      names(states) <- names(amal.deps$groups[[i]])
      states <- sapply(states, function(x) suppressWarnings(x[!is.na(as.numeric(x))]),  simplify = F )

      # Set initial Q matrices #
      M <- lapply(states, function(x) initQ(x,1) )

      # Decide between 'automatic' amalgamation (SMM) of 'check' amalgamation (evaluates if SMM or ED models) #
      if(mode == "auto"){ amal.deps$M[[i]] <- do.call(amaSMM, M) }
      if(mode == "check"){

        test <- sum(unlist(sapply(state.data[amal.deps$groups[[i]]], function(x) as.numeric(grep(x, pattern = "absent")) )))


        # Decide between SMM or ED models #
        if(test > 0){

          # ED-quality model #
          amal.deps$M[[i]] <- do.call(amaED, setNames(c(M, "ql"), c("Qc", rep("Qd", (length(M) - 1)), "type")))

          # OBS.: Should include ED-ap models (but complex, requires +3 traits) #

        }else{ amal.deps$M[[i]] <- do.call(amaSMM, M) } # SMM-ind models #

      }

    }else{

      # Single traits (non-dependent) #
      states <- unique(states)[order(unique(states))]
      states <- suppressWarnings(states[!is.na(as.numeric(states))])

      amal.deps$M[[i]] <- initQ(states,1)

    }

  }

  # Set names of Q matrices #
  names(amal.deps$M) <- names(amal.deps$groups)

  # Adjust names #
  names(amal.deps$state.data) <- names(amal.deps$M)

  # cleaning #
  amal.deps$groups <- amal.deps$groups[!amal.deps$drop]
  amal.deps$M <- amal.deps$M[!amal.deps$drop]
  amal.deps$state.data <- amal.deps$state.data[!amal.deps$drop]

  # Get states #
  amal.deps$states <- lapply(amal.deps$M, function(x) colnames(x) )

  # Return results #
  return(amal.deps)

}


#########################
##### OLD FUNCTIONS #####
#########################

#' Recodes a treedata object based on amalgamated characters.
#'
#' @import utils
#'
#' @title Recode a treedata object
#'
#' @description Recodes a treedata object based on amalgamated characters.
#'
#' @param td The treedata object to recode
#' @param traits The states that will be recoded
#' @param states The new character states for recoding
#' @param depstates Dependent states
#'
#' @examples
#' \dontrun{
#'     recode_td(treedata, traits = traits, states = states)
#' }
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
#' @title Recode Traits in a treedata Object
#'
#' @param td A 'treeplyr' treedata object with characters to recode
#' @param M A list produced by the function `amalgamate_deps`
#'
#' @examples
#' \dontrun{
#'     recode_traits(treedata, matrix)
#' }
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
#'
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
        ivs <- igraph::V(subgraph)
        dep_state <- get_dep_state(subgraph, ivs[ivs$name == node_name], bin_mats[[ancestor_v]], node_list) #trait_order)

        M[[i]] <- comb2matrices(bin_mats[[ancestor_v]], bin_mats[[node_name]], dependent.state = dep_state)
      }
      else { # otherwise, combine with the resulting matrix from the combination we just did
        # set dependent state
        ivs <- igraph::V(subgraph)
        dep_state <- get_dep_state(subgraph, ivs[ivs$name == node_name], M[[i - 1]], node_list) #trait_order)
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
#' Function for removing redundant dependency edges from our graph.
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


## Not exported (remove?) ##
empty.df = function(colnames)
{
    setNames(data.frame(matrix(ncol=length(colnames), nrow=0)),colnames)
}

                                        # this doesn't really work with lists...
addrow = function (df,row)
{
    colnames = names(df)
                                        #    print("addrow:")
                                        #    print(row)
                                        #    print(colnames)
    dfrow = data.frame(matrix(ncol=length(row),nrow=1,row), stringsAsFactors=FALSE)
    dfrow = setNames(dfrow,colnames)
    rbind(df,dfrow)
}

#
### Not exported (remove?) ##
##'
##'
#get.dependencies = function(uberon.filename, uterms, indirect=FALSE)
#{
#                                        # 1. Load uberon propagating different relationships
#    print("Loading ontology 'uberon.obo'")
#    uberon_is_a= get_ontology(uberon.filename,
#                              propagate_relationships = c("is_a"),
#                              extract_tags = 'minimal')
#
#    uberon_part_of = get_ontology(uberon.filename,
#                                  propagate_relationships = c("part_of"),
#                                  extract_tags = 'minimal')
#
#    uberon_develops_from = get_ontology(uberon.filename,
#                                        propagate_relationships = c("develops_from"),
#                                        extract_tags = 'minimal')
#
#    uberon = get_ontology(uberon.filename,
#                          propagate_relationships = c("is_a", "part_of","develops_from"),
#                          extract_tags = 'minimal')
#
#                                        # 2. Get names for the uberon terms
#    unames = c()
#    for(i in 1:length(uterms))
#    {
#        unames[i] = get_term_property(uberon,"name",uterms[i])
#    }
#
#                                        # 3. for each pair of terms, check dependencies
#    colnames = c("chr.id","dependent.chr.id","subrels")
#    df = empty.df(colnames)
#
#    if (indirect)
#    {
#        maybe.remove.indirect = function(x) {x}
#    }
#    else
#    {
#        maybe.remove.indirect = remove.indirect
#    }
#
#    D = maybe.remove.indirect(get_term_descendancy_matrix(uberon, uterms));
#    Disa = maybe.remove.indirect(get_term_descendancy_matrix(uberon_is_a, uterms));
#    Dpartof = maybe.remove.indirect(get_term_descendancy_matrix(uberon_part_of, uterms));
#    Ddevelopsfrom = maybe.remove.indirect(get_term_descendancy_matrix(uberon_develops_from, uterms));
#
#    for(i in 1:length(uterms))
#    {
#        for(j in 1:length(uterms))
#        {
#            if (i == j) next;
#
#            if (D[i,j])
#            {
#                cat(sprintf("%s %s\n",uterms[i],uterms[j]))
#
#                labels = c()
#
#                if (Disa[i,j])
#                {
#                    labels = c(labels,"is_a")
#                }
#                if (Dpartof[i,j])
#                {
#                    labels = c(labels, "part_of")
#                }
#                if (Ddevelopsfrom[i,j])
#                {
#                    labels = c(labels, "develops_from")
#                }
#
#                label = ""
#                if (length(labels) > 0)
#                {
#                    label = paste(labels,sep=":")
#                }
#
#                cat(sprintf("'%s' -> '%s'",unames[i],unames[j]),"\n")
#
#                # chr.id | chr.ancestor | label
#                df = addrow(df,c(uterms[i], uterms[j], label))
#
#                                        #                cat(sprintf("%s,0,%s,1\n",unames[j],unames[i]), file=outfile)
#                                        #                cat(sprintf("%s,1,%s,1\n",unames[j],unames[i]), file=outfile)
#
#
#                                        #                cat(edge, " ", attributes, "\n", file=dotfile)
#            }
#                                        # OK, so does i depend on j?
#        }
#    }
#
#    tuple(df,uterms,hashmap(uterms,unames),uberon)
#}


#############################################
##### Structured Markov Model functions #####
#############################################


# Tarasov et al. (2019): PARAMO #

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
#' g <- igraph::make_graph(c("A", "B", "B", "C", "C", "D"), directed = TRUE)
#' init_binary_matrix(g)
#' @importFrom igraph vcount
#' @importFrom igraph V
#' @importFrom igraph graph.edgelist
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
#' @importFrom igraph make_graph
#' @importFrom igraph topo_sort
#' @importFrom igraph ego
#' @examples
#' g <- igraph::make_graph(c("A", "B", "B", "C", "C", "D"), directed = TRUE)
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
#' @description Likely as independently evolving
#' @param list.matrices List of matrices
#' @return Matrix
#' @export
combNmatrices<-function(list.matrices){
  comb.matrix<-list.matrices[[1]]

  for (i in 1:(length(list.matrices)-1)){
    comb.matrix=comb2matrices(comb.matrix, list.matrices[[i+1]], dependent.state=NULL, diag.as = 0)
  }

  return(comb.matrix)
}


# Tarasov et al. (2019): PARAMO #

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


#' A wrapper for corHMM that fits Structured Hidden State Models.
#'
#' @param td A treeplyr treedata object
#' @param amalgamations A list produced by `amalgamate_deps`
#' @param ... Additional arguments passed to `corHMM::rayDISC`
#'
#' @details This function fits the models contained in `amalgamations$M`
#' to each of the characters in `td` using the function `corHMM::rayDISC`.
#'
#' @importFrom phytools bind.tip
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
    fits[[i]] <- corHMM::rayDISC(.phy, .dat, rate.mat = .M, model = "ARD", ...)
  }
  attributes(fits)$td <- td
  names(fits) <- colnames(td$dat)
  return(fits)
}


#' A wrapper for corHMM that draws stochastic character maps given a corHMM fit.
#'
#' @param fits A list produced by `amalgamated_fits_corHMM`
#' @param ... Additional arguments passed to `corHMM::makeSimmap`
#'
#' @details This function takes the models produced from `amalgamated_fits_corHMM` and
#' uses them to generate stochastic character maps.
#'
#' @importFrom phytools bind.tip
#' @importFrom phytools drop.tip.simmap
#'
#' @return A list of stochastic character maps for each character in the dataset
#'
#' @export
amalgamated_simmaps_corHMM <- function(fits, ...){
  td <- attributes(fits)$td
  trees <- list()
  for(i in 1:ncol(td$dat)){
    xx <- as.data.frame(cbind(td$phy$tip.label, as.character(td$dat[[i]])))
    Q <- fits[[i]]$solution
	phy <- phytools::bind.tip(td$phy, "...delete...", edge.length = 0, where=length(td$phy$tip.label)+1)
	new.row <- data.frame("...delete...", paste(0:(ncol(Q)-1), collapse="&"))
	colnames(new.row) <- colnames(xx)
	xx <- rbind(xx, new.row)
    stmp <- corHMM::makeSimmap(td$phy, xx, model=Q, rate.cat = 1, ...)

	trees[[i]] <- lapply(stmp, function(x) x <- phytools::drop.tip.simmap(x, "...delete...") )

  }
  return(trees)
}


#####################################
##### UPDATE: Sergei (2021)     #####
#####################################

## AMALGAM_CHAR ##


#' @name amaSMM
#' @aliases amaSMM_2Q
#' @rdname amaSMM
#'
#' @title SMM amalgamation of rate matrices
#' @description Amalgamating two (\code{amaSMM_2Q}) or several (\code{amaSMM}) rate matrices as  as independently (SMM-ind) or dependently (SMM-sw) evolving (Tarasov 2019, 2020). The parametric scheme of matrices is defined by integers, different integers
#' indicate different rate parameters
#' @param ... rate matrices
#' @param Qc rate matrix (controlling character if dependency is true)
#' @param Qd rate matrix (dependent character); if dependency is true then Qd depends on those states of Qc specified in controlling.state
#' @param controlling.state state(s) of Qc that switches on/off Qd matrix
#' @param name.sep separator for the state names
#' @param diag.as populates main diagonal elements
#' @param non.rate.as changes negative elements in the returned matrix to the specified value
#' @return matrix
#' @author Sergei Tarasov
#' @export
#'
#' @references Tarasov, S., 2019. Integration of anatomy ontologies and evo-devo using structured Markov models suggests a new framework for modeling discrete phenotypic traits. Systematic biology, 68(5), pp.698-716.
#' (\href{https://academic.oup.com/sysbio/article/68/5/698/5298740?login=true}{Read})
#' @references Tarasov, S., 2020. The invariant nature of a morphological character and character state: insights from gene regulatory networks. Systematic biology, 69(2), pp.392-400.
#' (\href{https://academic.oup.com/sysbio/article/69/2/392/5541792?login=true}{Read})
#'
#' @examples
#'C1 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
#'C2 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("o", "g"), c("o", "g")) )
#'C3 <-matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("w", "m"), c("w", "m")) )
#'amaSMM(C1, C2, C3)
amaSMM <- function(..., controlling.state = NULL, name.sep = "", diag.as = NULL, non.rate.as = NULL) {

  this_amaSMM_2Q <- function(Qc, Qd){

    amaSMM_2Q(Qc, Qd, controlling.state = controlling.state, name.sep = name.sep, diag.as = diag.as, non.rate.as = non.rate.as)

  }

  out <- Reduce(this_amaSMM_2Q, list(...))
  return(out)

}


#' @rdname amaSMM
#'
#' @export
#' @examples
#'Tl <- matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("a", "p"), c("a", "p")) )
#'C <-matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
#'amaSMM_2Q(Tl, C, controlling.state = NULL)
#'amaSMM_2Q(Tl, C, controlling.state = 2)
#'amaSMM(Tl, C, controlling.state = NULL)
amaSMM_2Q <- function(Qc, Qd, controlling.state = NULL, name.sep = "", diag.as = NULL, non.rate.as = NULL) {
  if (!is.null(controlling.state)) {
    matrix.diag <- rep(0, ncol(Qc))
    matrix.diag[controlling.state] <- 1
    matrix.diag <- diag(matrix.diag)
  }

  if (is.null(controlling.state)) {
    matrix.diag <- diag(nrow(Qc))
  }

  M_kr <- (Qc %x% diag(nrow(Qd))) + (matrix.diag %x% Qd)


  # getting colnames

  col <- paste(colnames(kronecker(Qc, diag(nrow(Qd)), make.dimnames = TRUE)),
               colnames(kronecker(diag(nrow(Qc)), Qd, make.dimnames = TRUE)),
               sep = ""
  )
  col <- gsub("::", name.sep, col, fixed = TRUE)

  # merging two names
  rownames(M_kr) <- colnames(M_kr) <- col
  if (!is.null(diag.as)) diag(M_kr) <- diag.as
  if (!is.null(non.rate.as)) M_kr[which(M_kr <= 0)] <- non.rate.as

  return(M_kr)
}



#' @name amaED
#' @aliases amaED_phi
#' @rdname amaED
#'
#' @title ED amalgamation of rate matrices
#' @description Amalgamating two matrices using embedded dependency (ED); \code{amaED} is the general function
#' that wraps \code{amaED_phi}. The parametric scheme of matrices is defined by integers, different integers
#' indicate different rate parameters
#'
#' @param Qc controlling character
#' @param Qd dependent character
#' @param type if \code{phi=NULL} specifies qualitative amalgamation (\code{ql}) or absent/present one (\code{ap}); for the latter absent state should be marked
#' with "*" (see the examples)
#' @param phi initial vector of the states of Qd when transiting from "absent" to "present" in Qc; \code{length(phi)} should be equal to \code{nrow(Qd)}
#' @param name.sep separator for the state names
#' @param ... other parameters for \code{amaED_phi}
#'
#' @details The ED amalgamation can be done automatically (\code{phi=NULL} and specified \code{type}) or manually if \code{phi} is provided.
#' For the automatic amalgamation \code{type='ap'}, the absent state in Qd should be marked with \code{'*'} as e.g., \code{'A*'}; in the automatic amalgamation  \code{phi} is
#' determined using \code{get_phi()}. When the state "present" in Qc controls more than one character (a/p and qualitative),
#' the qualitative character should be the second in amalgamation, e.g., \code{amaSMM(Qap, Qqual)}
#' @return matrix
#' @export
#'
#' @examples
#' # tail: absent (T*), present (T)
#' Tl <- matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("T*", "T"), c("T*", "T")) )
#' # color: red, blue
#' C <-matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
#' amaED(Tl, C, type=c("ql"), phi=NULL)
#' #
#' # tail armor
#' A <- matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("A*", "A"), c("A*", "A")) )
#' amaED(Tl, A, type=c("ap"), phi=NULL)
#' amaED(Tl, A, type=c("ap"), phi=c(1,0))
amaED <- function(Qc, Qd, type=c("ql", "ap"), phi=NULL, name.sep='', ...){
  phi1 <- phi
  type <- match.arg(type)

  if (is.null(phi1)){

    if (type=="ql"){
      cat("Using automatic qualitative type of intial vector phi.\n")
      phi.ql <- rep(1, nrow(Qd))
      out <- amaED_phi(Qc, Qd, phi=phi.ql, ...)
      return(out)
    }
    if (type=="ap"){
      cat("Using automatic a/p type of intial vector phi.\n")
      phi.ap <- get_phi(colnames(Qd), state.sep=name.sep)
      out <- amaED_phi(Qc, Qd, phi=phi.ap, ...)
      return(out)
    }

  } else {
    cat("Using the custom intial vector phi.\n")
    out <- amaED_phi(Qc, Qd, phi=phi1, ...)
    return(out)
  }

}


#' @rdname amaED
#'
#' @param non.rate.as changes negative elements in the returned matrix to the specified value
#' @param diag.as sets the main diagonal elements in the returned matrix
#' @export
#' @examples
#' Tl <- matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("T*", "T"), c("T*", "T")) )
#' C <-matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
#' amaED_phi(Tl, C, phi=c(1,1))
amaED_phi <- function(Qc, Qd, phi, name.sep = "", diag.as = NULL, non.rate.as = NULL){

  ap.label=rownames(Qc)

  Qnew <- Qd

  Qnew <- cbind(NA, Qnew)

  #-- check phi
  if (!is.numeric(phi)) stop("The intial vector phi should be numeric! Currently, phi=c(", paste0(phi, collapse=','), ")")
  if (length(phi)!=nrow(Qd)) stop("The length of intial vector phi should match nrow of Qd! Currently, length phi is ", length(phi))

  #--- report phi
  cat("The intial vector phi is phi=c(", paste0(phi, collapse=','), ")\n\n", sep='')
  #--


  phiQ <- c(0, Qc[1,2]*phi)
  Qnew <- rbind(phiQ, Qnew)

  Qnew[,1] <- c(rep(Qc[2,1], nrow(Qnew) ))
  diag(Qnew) <- 0
  diag(Qnew) <- -rowSums(Qnew)

  noms <- c(ap.label[1], paste0(ap.label[2], name.sep, colnames(Qd)) )
  rownames(Qnew) <- colnames(Qnew) <- noms

  #--
  if (!is.null(diag.as)) diag(Qnew) <- diag.as
  if (!is.null(non.rate.as)) Qnew[which(Qnew <= 0)] <- non.rate.as

  return(Qnew)

}


#' @title Determine initial vector phi for ED amalgamation
#' @description Given state names from rate matrix Qd constructs phi for ED amalgamation.
#' The absent state in Qd should be marked with \code{'*'} as e.g., \code{'A*'}
#'
#' @param mtnames state names from rate matrix
#' @param state.sep separator in amalgamated states
#' @return vector
#' @export
#'
#' @examples
#'A <- matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("A*", "A"), c("A*", "A")) )
#'get_phi(colnames(A))
get_phi <- function(mtnames, state.sep=''){
  nm <- strsplit(mtnames, state.sep)

  #n.stars <- sum(unlist(nm)=='*')
  n.stars <-lapply(nm, function(x) sum(x=='*'))
  n.stars <-unlist(n.stars)
  n.stars <-max(n.stars)

  if (n.stars==1){
    bool <- lapply(nm, function(x) ifelse(length(x)==1, FALSE, (x[2]=='*')) )
  } else if (n.stars>1) {
    bool <- lapply(nm, function(x) (x[2]=='*' & x[4]=='*') )
  } else {
    stop("Cannot construct phi automatically. Try manual construction of phi.")
  }

  #bool <- lapply(nm, function(x) ifelse(length(x)==1, FALSE, (x[2]=='*' & x[4]=='*') | (x[2]=='*')) )
  bool <- unlist(bool)

  if (all(bool==FALSE) | any(is.na(bool))) warning('No marks as * for absence-presence states are found in the Qd names.
  The automatic construction of the intial vector phi migh be incorrect. Consider specifying phi manually.')
  bool*1
}



#' @title Initialize custom rate matrix Q
#' @param char.state vector character states
#' @param rate.param vector of the rate parameters
#' @param diag.as sets the main diagonal elements; can be "negative", NA, or some value; "negative" returns negative row sum
#' @return matrix
#' @export
#' @examples
#' Q <- initQ(c('a', 'p'), c(1,2), diag.as = NA)
initQ <- function(char.state, rate.param, diag.as = "negative") {
  n.state <- length(char.state)
  Q <- matrix(ncol = n.state, nrow = n.state, byrow = TRUE)
  Q[xor(lower.tri(Q, diag = FALSE), upper.tri(Q, diag = FALSE))] <- rate.param
  Q <- t(Q)

  if (is.na(diag.as)){
    diag(Q) <- diag.as
  } else if (diag.as == "negative"){
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
  } else {
    diag(Q) <- diag.as
  }

  rownames(Q) <- colnames(Q) <- as.character(char.state)
  return(Q)
}




#' @title Converts rate matrix Q for using in RevBayes
#' @description Produces rate matrix that can be copy/paste into RevBayes script. The rate matrix can be represented as integers or symbols.
#'
#' @param Q rate matrix
#' @param round rate rounding
#' @param symb symbol for rates
#'
#' @return text string
#' @export
#'
#' @examples
#'Tl <- matrix(c(-1, 1, 1, -1), 2, 2, byrow = TRUE, dimnames =list( c("T*", "T"), c("T*", "T")) )
#'C <-matrix(c(-2, 2, 2, -2), 2, 2, byrow = TRUE, dimnames =list( c("r", "b"), c("r", "b")) )
#'Q <- amaSMM(Tl, C)
#'#
#'rb1 <- as_matrixRB(Q)
#'cat(rb1)
#'rb2 <- as_matrixRB(Q, symb='q')
#'cat(rb2)
as_matrixRB <- function(Q, round=2, symb=NULL){

  if (!is.null(symb)){
    out <- as_matrixRB_symbolic(Q, symb=symb)
    return(out)
  } else {
    out <- as_matrixRB_numeric(Q, round=round)
    return(out)
  }
}


# @title Converts rate matrix Q for using in RevBayes (entities as numbers)
as_matrixRB_numeric <- function(Q, round=2){
  diag(Q) <- 0
  #i=1
  out <- c()
  for(i in 1:nrow(Q)){
    x=Q[i,]
    x=format(round(x, round), nsmall = round)
    x=paste(x, collapse=', ')
    x=paste0('[', x, ']')
    out <- append(out, x)
  }

  out <- paste(out, collapse = ',\n ')
  out <- paste0('[', out, ']\n')
  return(out)

}


# @title Converts rate matrix Q for using in RevBayes (entities as symbols)
as_matrixRB_symbolic <- function(Q, symb='q'){
  round <- 0

  diag(Q) <- 0

  Q <- format(round(Q, round), nsmall = round)
  Qv <- c(Q)
  Qv <- sapply(Qv, function(x) ifelse(x=='0', '0.0', paste(symb, '[', x, ']',  sep='')) )
  Qnew <- matrix(Qv, nrow = nrow(Q), ncol = ncol(Q))

  out <- apply(Qnew, 1, function(x) paste0('[', paste0(x, collapse = ', '), ']'))

  out <- paste(out, collapse = ',\n ')
  out <- paste0('[', out, ']\n')
  #cat(out)
  return(out)

}



#####################################
##### UPDATE: Sergei (2022)     #####
#####################################

## HIDDEN_EXTENSION ##

#' @importFrom dplyr %>%
#' @importFrom MASS fractions
#' @importFrom numbers mLCM

#' @title Equal Rate Hidden Extension of rate matrix
#' @description Gives minimal equal rate hidden extension by finding least common multiple of rates.
#'
#' @param Q rate matrix
#' @return matrix
#' @export

#'
#' @examples
#'Q <- initQ(c(1, 2), c(.3,.2))
#'Qehe <- EHEtransform(Q)
#'
#'part_scheme=list(c(1, 2), c(3,4,5))
#'is_slumpable(Qehe, part_scheme)
EHEtransform<-function(Q)
{
  diag(Q)<-rep(0, nrow(Q))
  Q.vec<-c(Q)
  Q.vec.frac<-MASS::fractions(Q.vec, cycles = 10, max.denominator = 2000)
  fracs<-attr(Q.vec.frac[Q.vec.frac!=0],"fracs")
  # sometimes numbers without fraction appear like 1. add fractions
  fracs[fracs=="1"]<-"1/1"
  ####

  fracs<-strsplit(fracs, "/")
  vec.denom <- lapply(fracs, function (x) x[2]) %>% unlist() %>% as.numeric()

  lcm<-numbers::mLCM(vec.denom)
  lambda2<-1/lcm

  vec.num <- lapply(fracs, function (x) x[1]) %>% unlist() %>% as.numeric()
  theta.vec<-vec.num*(lcm/vec.denom)
  Theta<-Q
  Theta[Theta!=0]<-theta.vec

  # Expansion
  expn.by<-apply(Theta, 2, max)
  # entire row  can be 0 in theta, so change to 1
  expn.by[expn.by==0]<-1

  QD<-matrix(0, nrow = sum(expn.by), ncol = sum(expn.by), byrow = T)

  # make names
  names<-c()
  #i=2
  for (i in seq_along(expn.by))
  {
    names<-c(names, paste(i, ".", c(1:expn.by[i]), sep="") )
  }
  colnames(QD)<-rownames(QD)<-names

  # fill in lambda2
  # cumulative sum to find start of rows faster
  expn.by.cum<-cumsum(expn.by)
  expn.by.start<-c(1, expn.by.cum[-length(expn.by.cum)]+1 )

  for (i in 1:nrow(Theta))
  {
    for (j in 1:nrow(Theta))
    {
      if ( Theta[i, j]>0 )
      {
        QD[expn.by.start[i]:expn.by.cum[i],
           expn.by.start[j]:(expn.by.start[j]+(Theta[i, j]-1) )  ] <- lambda2
      }
    }
  }

  diag(QD)<-apply(QD, 1, function(x) sum(x))*-1
  return (QD)

}



#' @title Correlated Hyperspace Extension of rate matrix
#' @description Gives minimal Correlated Hyperspace Extension by first using EHEtransform() and then finding hypercube embedding.
#' So far, it works only for two-state rate matrices. Since there are many minimal CHE transformations for a given Q, it randomly generates one.
#'
#' @param Q rate matrix
#' @param threshold rounds the values in its first argument to the specified number of decimal places. It is needed for rate comparing rates using round().
#' @return a list that contains matrix and state names generated by latent characters

#' @export
#' @examples
#'Q <- initQ(c(1, 2), c(.3,.2))
#'tr <-CHEtransform(Q)
#'
#'Q.che <-tr$Q
#'part_scheme=list(c(1:4), c(5:8))
#'is_slumpable(Q.che, part_scheme)
CHEtransform <- function(Q, threshold=10){
  Qinit <- Q
  Q <- EHEtransform(Q)

  # get max number of hidden states per observable one
  cols <- colnames(Q)
  spl <- strsplit(cols, '\\.')
  obs.states <- lapply(spl, function(x) x[1]) %>% unlist %>% as.numeric()
  hid.states <- lapply(spl, function(x) x[2]) %>% unlist %>% as.numeric()
  n.obs <- unique(obs.states)
  hid.per.obs <- c()
  #i=1
  for (i in 1:max(n.obs)){
    hid.per.obs[i] <- max(hid.states[obs.states==i])
  }
  max.hid <- max(hid.per.obs)

  # calculate hyperspace
  qsin <- Q[Q>0][1]
  Q2 <- initQ(c(0:1), qsin)
  Q2.list <- vector(mode='list', length = max.hid)
  Q2.list <-lapply(Q2.list, function(x) Q2)

  #amaSMM(rep(Q2,2))
  Qche <- Reduce(amaSMM, Q2.list)
  cols.che <- colnames(Qche)
  # calculates sum to divide states into odd and even
  slist <- strsplit(cols.che, '')
  ssum <- lapply(slist, function(x) as.numeric(x) %>% sum)%>% unlist
  # color #1, 0 & even
  c1 <- which((ssum%%2)==0)
  # color #2, odd
  c2 <- which((ssum%%2)!=0)
  v <- c(c1,c2)
  Qche <- Qche[v,v]

  # make Qche correlated
  nstates.block <- 2^max.hid/2

  # block 12
  block12 <- c((nstates.block+1):ncol(Qche))
  block11 <- c(1:nstates.block)
  Qbl <- Qche[block11, block12]
  rate.focal <- Qinit[1,2]

  # using threshold since sometimes numbers do not match
  bool <- round(rowSums(Qbl), threshold)==rate.focal
  if (all(bool)==FALSE) {
    Qbl[Qbl==0] <- NA
    Qbl[!is.na(Qbl)] <- 0
    n.cells <- rate.focal/qsin

    # randomly populate cells
    Qbl <- apply(Qbl, 1, function(x) { smp <- sample(c(1:max.hid), n.cells);
    x[!is.na(x)][smp] <- qsin; x } ) %>% t()

    Qbl[is.na(Qbl)] <- 0
    Qche[block11, block12] <- Qbl
  }

  # block 21
  Qbl <- Qche[block12, block11]
  rate.focal <- Qinit[2,1]

  # using threshold since sometimes numbers do not match
  bool <- round(rowSums(Qbl), threshold)==rate.focal
  if (all(bool)==FALSE) {
    Qbl[Qbl==0] <- NA
    Qbl[!is.na(Qbl)] <- 0
    n.cells <- rate.focal/qsin

    # randomly populate cells
    Qbl <- apply(Qbl, 1, function(x) { smp <- sample(c(1:max.hid), n.cells);
    x[!is.na(x)][smp] <- qsin; x } ) %>% t()

    Qbl[is.na(Qbl)] <- 0
    Qche[block12, block11] <- Qbl
  }

  n1 <- paste0(rep(1, nstates.block), '.', 1:nstates.block)
  n2 <- paste0(rep(2, nstates.block), '.', 1:nstates.block)

  old.names <- colnames(Qche)
  colnames(Qche) <- row.names(Qche) <- c(n1,n2)
  diag(Qche) <- 0
  diag(Qche) <- -rowSums(Qche)

  return(list(Q=Qche, corr.states=old.names))
}



#' @title Convert an instance of rate matrix to model
#' @param Q an instance of rate matrix with cells having specific rate values
#' @param diag.as sets the main diagonal elements; can be "negative", NA, or some value; "negative" returns negative row sum
#' @return matrix where numbers correspond to different rate parameters
#' @export

#' @examples
#'Q <- initQ(c(1, 2), c(.3,.2))
#'Q2model(Q)
Q2model<-function(Q, diag.as=NA)
{
  uniq <- Q[Q>0] %>% unique()
  for (i in 1:length(uniq))
  {
    Q[Q==uniq[i]]<-i
  }

  Q[Q<=0]<-NA

  if (is.na(diag.as)){
    diag(Q) <- NA
  } else if (diag.as == "negative"){
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
  } else {
    diag(Q) <- diag.as
  }

  return (Q)
}


##################
#
# check strong lumpability
#
###########

#' @title Tests if a rate matrix is strongly lumpable
#' @description Given rate matrix Q and state partition, checks if Q is strongly lumpable
#' @param Q a rate matrix with cells having specific rate values
#' @param part_scheme partition of state in Q represented as a list
#' @return Boolean value
#' @export

#' @examples
#'Q <- initQ(c(1, 2), c(.3,.2))
#'Q.h <- EHEtransform(Q)
#'
#'part_scheme=list(c(1, 2), c(3,4,5))
#'is_slumpable(Q.h, part_scheme)
#'
#'part_scheme=list(1, c(2, 3, 4, 5))
#'is_slumpable(Q.h, part_scheme)
#'
#'Q.h <-CHEtransform(Q)$Q
#'
#'part_scheme=list(c(1:4), c(5:8))
#'is_slumpable(Q.h, part_scheme)
#'
#'part_scheme=list(c(1:3), c(4:8))
#'is_slumpable(Q.h, part_scheme)
is_slumpable<-function(Q, part_scheme){

  # check if num of chars in scheme match that of Q
  if (length(unlist(part_scheme))!=ncol(Q))
    stop("The number of states in part_scheme and Q does not match")

  # normilize diag
  diag(Q)<-apply(Q, 1, function(x) sum(x))*-1

  Nper.part=lapply(part_scheme, length)%>%unlist
  stat2M=matrix(0, nrow=length(Nper.part), ncol=ncol(Q))
  for (i in 1:nrow(stat2M))
    stat2M[i, part_scheme[[i]]]<-1

  M_rows=Q %*% t(stat2M)
  tru.vals=c()

  for (i in 1:length(part_scheme)){
    for (j in 1:ncol(M_rows)){
      tru.vals=c(tru.vals, length(unique(M_rows[part_scheme[[i]],j]))<2)
    }
  }
  # if this is false then matrix is not lumpable
  return(all(tru.vals))
}
