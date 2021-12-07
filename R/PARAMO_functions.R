

##############################
##### RevBayes functions #####
##############################


#' Data object for Rev scripts.
#'
#' @export
ratemat1 <- function() {
  '\nfor (i in 1:NUM_STATES) {\nfor (j in 1:NUM_STATES) {\nrates[i][j] <-0.0\n}\n}\n#rate prior\nr1 ~ dnExp(20)\nr2 ~ dnExp(20)\n\nmoves[++mvi] = mvScale(r1, lambda=1, tune=true, weight=2)\nmoves[++mvi] = mvScale(r2, lambda=1, tune=true, weight=2)
\n\n# place rate categories into matrix\nrates[2][1]:=r1\nrates[1][2]:=r2\n\n\nrate_matrix := fnFreeK(transition_rates=rates, rescaled=false, matrixExponentialMethod="eigen")\n\nroot_freq <- simplex(1, 1)\n\n'
  
}


#' Function to write Ontology CTMC models as a .Rev file for execution in RevBayes.
#'
#' @export
writeRevOntologyModels <- function(td, M, dir, dirW, dirR, dirD) {
  MT <- as.data.frame(td$dat)
  colnames(MT) <- gsub(" ", "_", colnames(MT))
  rownames(MT) <- td$phy$tip.label
  for (i in 1:ncol(MT)){
    C.rev<-MT[,i]
    C.rev<-gsub("&", " ", C.rev)
    o <- order(nchar(C.rev))
    
    out<-cbind(rownames(MT), C.rev)
    out <- out[o,]
    write.table(file=paste0(dirD, colnames(MT[i]), ".char"), out, quote=F, sep=" ", 
                row.names=F, col.names=F)
  }
  
  data(ParamoRevTemplate)
  
  
  for (i in 1:ncol(MT)){
    fl.in  <- ParamoRevTemplate
    fl.in  <- gsub(pattern = "Hymenoptera_br_resolved", replace = "fishtree",
                   x = fl.in)
    fl.in  <- gsub(pattern = "@analysis_name@", replace = paste0(colnames(MT[i])),
                   x = fl.in)
    
    fl.in <- gsub(pattern = "@chrs_2_read@", 
                  replace = paste0("data/", colnames(MT[i]), ".char"), x = fl.in)
    
    if(colnames(MT)[i] %in% gsub(" ", "_", names(M))){
      in.rev<-Mk_Rev(M[[gsub("_", " ", colnames(MT)[i])]])
      
      fl.in <- gsub(pattern = "@numstates@", 
                    replace = as.character(max(dim(M[[gsub("_", " ", colnames(MT)[i])]]))), x = fl.in)
      
      fl.in <- gsub(pattern = "@ratematrix@", 
                    replace = in.rev, x = fl.in)
    } else {
      
      fl.in <- gsub(pattern = "@numstates@", 
                    replace = "2", x = fl.in)
      
      fl.in <- gsub(pattern = "@ratematrix@", 
                    replace = ratemat1(), x = fl.in)
    }
    
    cat(file=paste0(dir, colnames(MT[i]), ".Rev"), sep="\n", fl.in)
  }
  
  
}


########################################
##### Stochastic Mapping functions #####
########################################


#' Make branch lengths and maps equal, use tree.tmp with rounded (3) edge.length.
#'
#' @export
make_tree_eq<-function(tree.tmp, target.tr, round=3){
  target.tr$edge.length<-tree.tmp$edge.length
  
  target.tr$maps<-lapply(target.tr$maps, function(x) round(x,3) )
  Maps.targ<- unlist(lapply(target.tr$maps, function(x) sum(x)))
  
  #k<-tree.tmp$edge.length-Maps.targ
  k<-tree.tmp$edge.length/Maps.targ
  k<-as.list(k)
  
  maps.out<-mapply(function(x,y) 
  {x*y },
  x=target.tr$maps, y=k )
  
  target.tr$maps<-lapply(maps.out, function(x) round(x, round))
  return(target.tr)
}


#' Reading unsummarized Stochastic Map files from RevBayes.
#'
#' @param file file
#' @param start start from tree
#' @param end end with tree
#' @param save save to file. if NULL reads in R
#' @export
read_Simmap_Rev<-function(file, start=1, end=1, save=NULL){
  
  skip=start+2
  max2read=end-start+1
  
  text <- scan(file=file, sep = "\n", what = "character", skip=skip, nlines=max2read)
  
  trees<-c()
  for (i in 1:length(text)){
    
    #trees[i]<-strsplit(text[i], "\\}\t\\(")[[1]][2]
    
    ss=regexpr("\\}\t\\(",  text[i])[1]
    trees[i]<-substring(text[i], first=ss+2)
  }
  
  if (is.null(save)){
    return(trees)
  }else{
    
    cat(trees, file=save, sep="\n")
    print(paste0("Tree(s) are saved to ", save))
  }
}


#' Reading unsummarized Stochastic Map files from RevBayes for one tree.
#'
#' @export
discr_Simmap<-function(tree, res){
  
  steps <- 0:res/res * max(phytools:::nodeHeights(tree))
  H <- phytools:::nodeHeights(tree)
  maps.n <- vector(mode = "list", length = nrow(tree$edge))
  
  # i=170
  for (i in 1:nrow(tree$edge)) {
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]), 
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) -  H[i, 1]
    
    
    TR<-cumsum(tree$maps[[i]])
    # TR[length(TR)]<-YY[nrow(YY), 2] # this to make the length equal as it sometiems does not hold
    # YY[,1]-YY[,2]

    #sprintf("%.54f", c(TR[length(TR)], YY[nrow(YY), 2]) )
    # TR[1]==TR[2]
    # all.equal(TR[1], TR[2])
    # all.equal(TR)
    # duplicated(TR)
    #length(int.out)
    #all.equal(TR[length(TR)], YY[nrow(YY), 2])
    #TR[length(TR)]==YY[nrow(YY), 2]
    
    #TR[length(TR)]==YY[nrow(YY), 2]
    int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE, all.inside = TRUE)
    #int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE)
    #findInterval(seq(0.1, 4, .1), c(0, 0.5, 0.7, 1.5, 1.6, 4 ), left.open=T, rightmost.closed = F)
    maps.n[[i]]<-setNames(YY[,2]-YY[,1], names(tree$maps[[i]])[int.out])
  }
  tree$maps<-maps.n
  return(tree)        
}


#' Reading unsummarized Stochastic Map files from RevBayes.
#'
#' @export
discr_Simmap_all<-function(tree, res){
  
  if (class(tree)[1]=="simmap") {
    tree<-discr_Simmap(tree, res)
  }
  
  if (class(tree)[1]=="multiSimmap") {
    
    for (j in 1:length(tree)){
      tree[[j]]<-discr_Simmap(tree[[j]], res)
    }
  }
  return(tree)
}


#' Function for processing RevBayes stochastic maps (V1).
#'
#' @export
prepareMapsRev <- function(td, discretization_level=100, start_tree=1, end_tree=2) {
  characters <- colnames(td$dat)
  characters <- gsub(" ", "_", characters)

  # Read a sample of 2 maps from .stm files and save them in the proper format .stmR
  
  for (i in 1:length(characters))
  {
    .tree<-read_Simmap_Rev(paste0(dirR, characters[i], ".stm"),
                           start=start_tree, end=end_tree,
                           save = NULL) 
    tree  <- phytools::read.simmap(text=.tree, format="phylip")
    
    
    phytools::write.simmap(tree, file=paste0(dirW, characters[i], ".stmR"))
  }
  
  # Read stmR, discretize maps, and save each map as a separate rds file; 
  # All rds files for a chracter are stored in a zip archive.

  
  for (i in 1:length(c))
  { 
    # read in undescritezed trees
    #print(paste0("Reading ", characters[i]))
    sim=phytools::read.simmap(file=paste0(dirW, characters[i], ".stmR"), format="phylip")
    
    # discretize trees by looping over sample and saving as rds
    
    for (j in 1:length(sim)){
      tryCatch({
        
        #print(paste0("Discretizing tree ", j))
        
        ## errors with na
        
        # make trees equal with template
        sim.d<-make_tree_eq(td$phy, sim[[j]], round=5)
        
        #sim.d<-discr_Simmap_all(sim[[j]], 1000)
        sim.d<-discr_Simmap_all(sim.d, discretization_level)
        
        saveRDS(sim.d, file =  paste0(dirW,characters[i], "_", j, ".rds") )
        
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        #errors<-rbind(errors, c(ii,jj))
      }  )
      
    } 
    
    # putting rds files into archive
    files<-paste0(dirW, characters[i], "_", c(1:length(sim)), ".rds")
    zip(paste0(dirW, characters[i], ".zip"), files=files)
    file.remove(files)
    
  }
  
  # close connections
  showConnections (all=T)
  closeAllConnections()
  
}


#' Function for processing RevBayes stochastic maps (V2).
#'
#' @export
prepareMapsRayDISC <- function(td, simmaps, discretization_level=100) {
  characters <- colnames(td$dat)
  characters <- gsub(" ", "_", characters)
  names(simmaps) <- gsub(" ", "_", names(simmaps))

  # Read a sample of 2 maps from .stm files and save them in the proper format .stmR

  #trees <- list()
  #for (i in 1:length(characters))
  #{
  #  tree[[i]]  <- simmaps$trees[i]
  #  write.simmap(tree, file=paste0(dirW, characters[i], ".stmR"))
  #}
  

  # Read stmR, discretize maps, and save each map as a separate rds file; 
  # All rds files for a chracter are stored in a zip archive

  
  if(!"multiSimmap" %in% class(simmaps[[1]])){
    DISC_trees <- list()
    for (i in 1:length(characters)){ 
      tryCatch({
        
        #print(paste0("Discretizing tree ", j))
        
        ## errors with na
        
        # make trees equal with template
        sim.d<-make_tree_eq(td$phy, simmaps[[characters[i]]], round=5)
        
        #sim.d<-discr_Simmap_all(sim[[j]], 1000)
        sim.d<-discr_Simmap_all(sim.d, discretization_level)
        DISC_trees[[i]] <- list(sim.d)
        #saveRDS(sim.d, file =  paste0(dirW,characters[i], "_", j, ".rds") )
        
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        #errors<-rbind(errors, c(ii,jj))
      }  )
    }
} else {
  DISC_trees <- lapply(1:length(characters), function(x) list())
  { 
    for (i in 1:length(characters)){ 
    tryCatch({
      
      #print(paste0("Discretizing tree ", j))
      
      ## errors with na

      for(j in 1:length(simmaps[[i]])){
        sim.d<-make_tree_eq(td$phy, simmaps[[i]][[j]], round=5)
        
        #sim.d<-discr_Simmap_all(sim[[j]], 1000)
        sim.d<-discr_Simmap_all(sim.d, discretization_level)
        DISC_trees[[i]][[j]] <- sim.d
      }
      # make trees equal with template
      #saveRDS(sim.d, file =  paste0(dirW,characters[i], "_", j, ".rds") )
    }, error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      #errors<-rbind(errors, c(ii,jj))
    }  )
}
  }
  }
  names(DISC_trees) <- characters

    # putting rds files into archive
    #files<-paste0(dirW, characters[i], "_", c(1:length(sim)), ".rds")
    #zip(paste0(dirW, characters[i], ".zip"), files=files)
    #file.remove(files)
    
  # close connections
  #showConnections (all=T)
  #closeAllConnections()
  return(DISC_trees)
}


## Not exported ##
stack_stm<-function(stm.list){
  M<-lapply(stm.list, function(x) x$maps)
  M<-lapply(M, function(x) lapply(x, function(y) names(y)))
  M<-Reduce(stack2, M)
  
  M.out<-mapply(function(x,y) 
  {setNames(x, y) },
  x=stm.list[[1]]$maps, y=M )
  
  out<-stm.list[[1]]
  out$maps<-M.out
  return(out)
  
}


## Not exported ##
# Stack two discrete stm's lists; x,y are the list of state names (i.e. maps).
stack2<-function(x,y){
  mapply(function(x,y) 
  {paste(x,y, sep="") },
  x=x, y=y )
}


#' Final stack of maps for a set of stochastic maps stored in a directory.
#'
#' @param cc chars id to stack
#' @param ntrees number of trees to stack
#' @param dirW directory for zip file
#' @return A list of stacked stochastic character maps
#' @export
paramo<-function(cc, ntrees=10, dirW=c("") )
{
  tr<-vector("list", ntrees)
  for (i in 1:ntrees){
    
    fl<-paste0(cc, "_", i, ".rds")  
    
    stack.L<-vector("list", length(fl))
    
    for (j in 1:length(fl)){
      
      print(paste0("Reading ", paste0(cc[j], ".zip"), " and ", fl[j]))
      con<-unz(paste0(dirW, cc[j], ".zip"), filename=paste0(dirW, fl[j]) )
      con2 <- gzcon(con)
      stack.L[[j]] <- readRDS(con2)
      close(con)
    }
    
    tr[[i]]<- stack_stm(stack.L)
  }
  return(tr)
}


#' Final stack of maps for a set of stochastic maps stored in a list.
#'
#' @param cc chars id to stack
#' @param tree.list Named list with stochastic character maps
#' @param ntrees number of trees to stack
#' @return A list of stacked stochastic character maps
#' @export
paramo.list<-function(cc, tree.list, ntrees=1)
{
  tr<-vector("list", ntrees)
  ncharacters <- length(cc)
  cc <- gsub(" ", "_", cc)
  for (i in 1:ntrees){
    stack.L<-vector("list", length(cc))
    for (j in 1:ncharacters){
      stack.L[[j]] <- tree.list[[cc[j]]][[i]]
      #names(stack.L[[j]]) <- cc[j]
    }
    tr[[i]]<- stack_stm(stack.L)
  }
  return(tr)
}


#####################################
##### General utility functions #####
#####################################


#' Function for aggregating characters under a specified set of terms (for example, body regions).
#'
#' @export
RAC_query <- function(char_info, ONT, names){
  c <-  char_info$ID
  c <- gsub(" ", "_", c)
  
  annot <- as.list(as.character(char_info$IRI))
  names(annot) <- as.character(char_info$ID)
  ONT$terms_selected_id <- annot
  
  
  parts <- do.call(rbind,lapply(names,  pk_get_iri, as="uberon"))
  parts$IRI <- sapply(parts[,1], strip_IRI)
  levelA <- setNames(parts$IRI, names)     
  
  res <- lapply(levelA, function(x)
    get_descendants_chars(ONT, annotations="manual", terms=x)  )
  
  cat("\nAggregations by :\n")
  print(res)
  return(res)
}


#' Function to convert a list of character annotations to matrix.
#'
#' @title Convert list to edge matrix
#' @description Takes list of charater annotations and creates an edge matrix comprising two columns: from and to.
#' The list to table conversion can be done using ldply function from plyr package: plyr::ldply(list, rbind).
#' @param annotated.char.list Character list with ontology annotations.
#' @param col_order_inverse The default creates the first columns consisting if character IDs and the second columns consisting of ontology annatotaions.
#' The inverse order changes the columns order.
#' @return Two-column matrix.
#' @examples
#' annot_list<-list(`CHAR:1`=c("HAO:0000933", "HAO:0000958"), `CHAR:2`=c("HAO:0000833", "HAO:0000258"))
#' list2edges(annot_list)
#' # attache plyr package and run
#' # ldply(annot_list, rbind)
#' @export
list2edges<-function(annotated.char.list, col_order_inverse=F){
  annotated.vec=setNames(unlist(annotated.char.list, use.names=F),rep(names(annotated.char.list), lengths(annotated.char.list)))
  if (col_order_inverse==T){
    edge.matrix=cbind(unname(annotated.vec), names(annotated.vec))
  } else
    edge.matrix=cbind(names(annotated.vec), unname(annotated.vec))
  return(edge.matrix)
}


#' Function to convert a table of character annotations to list. 
#'
#' @title Convert table to list
#' @description Takes table where each row consists of character number and ontology annotations and returns a list.
#' Each character is assigned its own ID CHAR:XXXX
#' @param table A character table with annotations.
#' @param id_col A column ID corresponding to character
#' @param descendants_cols IDs of columns corresponding to character annotations
#' @return The list.
#' @examples
#' # converting Sharkey_2011 data set to list of characater states
#' list_data<-table2list(Sharkey_2011)
#' @export
table2list<-function(table, id_col=c(1), descendants_cols=c(2:ncol(table))){
  annotated.char.list=list()
  for (i in 1:nrow(table)) {
    annotated.char.list[[i]]=c(na.omit(as.character(table[i,descendants_cols])))
  }
  #names(annotated.char.list)<-paste("CHAR:", table[,id_col], sep="")
  return(annotated.char.list)
}


#' Function to extract IDS from ontology terms.
#'
#' @title Get IDs for ontology names
#' @description Returns IDs of ontology terms given terms' names
#' @param vec_name names od terms
#' @param ontology ontology
#' @param names use element name
#' @return vector of IDs.
#' @examples
#' vec_name=c("ventral mesofurco-profurcal muscle", "anatomical entity")
#' get_onto_id(vec_name, HAO)
#' @export
get_onto_id<-function(vec_name, ontology, names=F){
  match_vec<-match(unlist(vec_name, use.names = FALSE), ontology$name)
  ids=names(ontology$name)[match_vec]
  if (names==T) {names(ids)<-ontology$name[match_vec]}
  return(ids)
}


#' Function to get all characters annotated with terms descending from a selected ontology term.
#
#' @title Get characters that descendants of selected ontology term
#' @description Returns all characters located (associated) with a given ontology terms
#' @param ontology ontology_index object.
#' @param annotations which annotations to use: "auto" means automatic annotations, "manual" means manual ones.
#' Alternatively, any othe list element containing annotations can be specified.
#' @param terms IDs of ontology terms for which descendants are queried.
#' @param ... other parameters for ontologyIndex::get_descendants() function
#' @return The vector of character IDs.
#' @examples
#' ontology<-HAO
#' ontology$terms_selected_id<-list(`CHAR:1`=c("HAO:0000653"), `CHAR:2`=c("HAO:0000653"))
#' get_descendants_chars(ontology, annotations="manual", "HAO:0000653")
#' @export
get_descendants_chars<-function(ontology, annotations="auto", terms, ...){
  
  if (is.list(annotations)){
    annot_list<-annotations # specify your annotation list
  } else {
    
    if (annotations=="auto"){
      annot_list<-ontology$auto_annot_characters
    }
    if (annotations=="manual"){
      annot_list<-ontology$terms_selected_id
    }
  }
  
  
  onto_chars_list=list2edges(annot_list, col_order_inverse=T)
  descen<-unique(onto_chars_list[,2][onto_chars_list[,1] %in%
                                       ontologyIndex::get_descendants(ontology=ontology, roots=terms, ...)])
  return(descen)
}