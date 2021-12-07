#' Function to display the number of matching species and average value across them from an Ontotrace matrix
#' @export
print_coverage <- function(x){
  coverage <- apply(x, 2, function(x) sum(!is.na(x)))
  average <- sapply(x, function(x) mean(as.numeric(na.omit(x[x %in% c("0", "1")])), na.rm=TRUE))
  cover <- cbind(coverage, average)
  tmp <- dplyr::filter(data.frame(traits=rownames(cover), cover), coverage > 0, average < 1, average > 0) %>% arrange(., desc(coverage))
  print(tmp)
}

#' Utility function for stripping off url from IRI
#' @export
strip_IRI <- function(x){
  x <- gsub("http://purl.obolibrary.org/obo/", "", x)
  x <- gsub("_", ":", x)
  return(x)
}

#' Utility function for cleaning up character data table after amalgamating characters
#' @export
dropDependentTraits <- function(char_info, dep.mat, td){
  char_info_comb <- char_info[which(apply(dep.mat, 1, sum, na.rm=TRUE)==0), c(1,5)]
  new.traits <- colnames(td$dat)
  old.traits <- sapply(new.traits, function(x) strsplit(x, "+", fixed=TRUE)[[1]][1])
  trait.trans <- setNames(new.traits, old.traits)
  char_info_comb$ID <- unname(trait.trans[as.character(char_info_comb$ID)])
  return(char_info_comb)
}

#' Utility function for filtering based on missing traits and taxa
#' @export
filter_coverage <- function(td, traits=0, taxa=0){
  taxa_coverage <- apply(td$dat, 1, function(x) mean(as.numeric(!is.na(x))))
  trait_coverage <- apply(td$dat, 2, function(x) mean(as.numeric(!is.na(x))))
  td <- dplyr::filter(td, taxa_coverage > taxa)
  td <- dplyr::select(td, which(trait_coverage > traits))
  return(td)
}

#' Data object for Rev scripts
#' @export
ratemat1 <- function() {
  '\nfor (i in 1:NUM_STATES) {\nfor (j in 1:NUM_STATES) {\nrates[i][j] <-0.0\n}\n}\n#rate prior\nr1 ~ dnExp(20)\nr2 ~ dnExp(20)\n\nmoves[++mvi] = mvScale(r1, lambda=1, tune=true, weight=2)\nmoves[++mvi] = mvScale(r2, lambda=1, tune=true, weight=2)
\n\n# place rate categories into matrix\nrates[2][1]:=r1\nrates[1][2]:=r2\n\n\nrate_matrix := fnFreeK(transition_rates=rates, rescaled=false, matrixExponentialMethod="eigen")\n\nroot_freq <- simplex(1, 1)\n\n'
  
}

#' recodes a treedata object based on amalgamated characters
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

#' Plots a heatmap along with a phylogeny and trait tree
#' @export
ontologyHeatMap <- function(td, njt, start=3, margs=c(0.2, 0.25), ...){
  #vals <- na.omit(unique(do.call(c, lapply(3:ncol(td$dat), function(x) unique(as.character(td$dat[[x]]))))))
  
  X <- do.call(cbind, lapply(start:ncol(td$dat), function(x) as.numeric(recode(td$dat[[x]], "0 and 1"=0.5, "1 and 0"=0.5, "1"=1, "0"=0, "2"=2, "3"=3))))
  colnames(X) <- colnames(td$dat)[start:ncol(td$dat)]
  
  .vals <- sort(na.omit(unique(as.vector(X))))
  dimx <- dim(X)
  
  tree1 <- phytools::force.ultrametric(td$phy,method = "extend")
  tree_ord <- attr(tree1, "order")
  
  if(!is.null(njt)){
    X <- X[,njt$tip.label[njt$edge[njt$edge[,2] <= length(njt$tip.label),2]]]
    tree2 <- njt
    tree2 <- ape::chronopl(njt, 1)
    tree2$edge.length <- tree2$edge.length/(max(branching.times(tree2)))*margs[2]*dimx[1]
    png(tempfile())
    invisible(h2 <- plot(tree2, plot = FALSE, direction = "downwards", show.tip.label=FALSE))
    dev.off()
  } else{
    h2 <- list(x.lim=c(1,dimx[2]+1), y.lim=c(0,0.2*dimx[1]))
  }
  
  #alters the length of the tips of the trees
  tree1$edge.length <- tree1$edge.length/(max(branching.times(tree1)))*margs[1]*dimx[2]
  
  #Changes the direction of the top plot
  png(tempfile())
  invisible(h1 <- plot(tree1, plot = FALSE, cex=0.5))
  dev.off()
  
  # adjustible color palette for the plot and legend
  colors <- c("#ffeaa7","#fab1a0", "#e17055")
  
  
  # this is all setting boundaries for the different plots and combining them into one image
  par(mar = c(0,0,0,0))
  plot(0,0, type = 'n', xlim = c(0,h1$x.lim[2]+h2$x.lim[2]), ylim=c(0,h1$y.lim[2]+h2$y.lim[2]))
  
  image(seq(h1$x.lim[2]+1,h1$x.lim[2]+h2$x.lim[2], length.out=ncol(X)), seq(1, h1$y.lim[2], length.out=nrow(X)), t(X), xlim=c(1+h1$x.lim[2],h1$x.lim[2]+h2$x.lim[2]+1) ,ylim=c(0, h1$y.lim[2]-1), add=TRUE, col=colors)
  
  legend(0, (h1$y.lim[2]+h2$y.lim[2])*.99, legend=.vals ,pch=22, pt.bg=colors)
  
  par(new = TRUE)
  
  plot(tree1, x.lim=c(0,(1+margs[2])*(h2$x.lim[2]+h1$x.lim[1])), y.lim=c(0,h1$y.lim[2]+h2$y.lim[2]),...)
  
  if(!is.null(njt)){
    par(new = TRUE)
    plot(tree2, direction = "downwards", x.lim=c(-h1$x.lim[2],h2$x.lim[2]), y.lim=c((-h1$y.lim[2])-0.01*dimx[1],h2$y.lim[2]))
  }
  
  return ()
}

#' Makes a trait tree using semantic similarity
#' @export
makeTraitTree <- function (td, skip=1:2){
  traits <- colnames(td$dat)
  traits <- traits[-(1:2)] #delete otu data
  traits
  
  semanticSimilarityMatrix <- jaccard_similarity(terms = traits, .colnames = "label", .labels = traits)
  
  #rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- traits
  
  #Neighbor-joining tree of SS matrix
  njt <- nj(1-semanticSimilarityMatrix)
  
  return (njt)
}

#' Function for processing revbayes stochastic maps
#' @export
prepareMapsRev <- function(td, discretization_level=100, start_tree=1, end_tree=2) {
  characters <- colnames(td$dat)
  characters <- gsub(" ", "_", characters)
  #####################################
  # Read a sample of 2 maps from .stm files and save them in the proper format .stmR
  #####################################
  
  for (i in 1:length(characters))
  {
    .tree<-read_Simmap_Rev(paste0(dirR, characters[i], ".stm"),
                           start=start_tree, end=end_tree,
                           save = NULL) 
    tree  <- phytools::read.simmap(text=.tree, format="phylip")
    
    
    phytools::write.simmap(tree, file=paste0(dirW, characters[i], ".stmR"))
  }
  ##########
  
  #####################################
  # Read stmR, discretize maps, and save each map as a separate rds file; 
  #all rds filea for a chracter are stored in a zip archive
  #####################################
  
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
        
        ##
        
        ##### make trees equal with template
        sim.d<-make_tree_eq(td$phy, sim[[j]], round=5)
        ###
        
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

#' Function for processing revbayes stochastic maps
#' @export
prepareMapsRayDISC <- function(td, simmaps, discretization_level=100) {
  characters <- colnames(td$dat)
  characters <- gsub(" ", "_", characters)
  names(simmaps) <- gsub(" ", "_", names(simmaps))
  #####################################
  # Read a sample of 2 maps from .stm files and save them in the proper format .stmR
  #####################################
  #trees <- list()
  #for (i in 1:length(characters))
  #{
  #  tree[[i]]  <- simmaps$trees[i]
  #  write.simmap(tree, file=paste0(dirW, characters[i], ".stmR"))
  #}
  ##########
  
  #####################################
  # Read stmR, discretize maps, and save each map as a separate rds file; 
  #all rds filea for a chracter are stored in a zip archive
  #####################################
  
  if(!"multiSimmap" %in% class(simmaps[[1]])){
    DISC_trees <- list()
    for (i in 1:length(characters)){ 
      tryCatch({
        
        #print(paste0("Discretizing tree ", j))
        
        ## errors with na
        ##
        
        ##### make trees equal with template
        sim.d<-make_tree_eq(td$phy, simmaps[[characters[i]]], round=5)
        ###
        
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
      ##
      for(j in 1:length(simmaps[[i]])){
        sim.d<-make_tree_eq(td$phy, simmaps[[i]][[j]], round=5)
        ###
        
        #sim.d<-discr_Simmap_all(sim[[j]], 1000)
        sim.d<-discr_Simmap_all(sim.d, discretization_level)
        DISC_trees[[i]][[j]] <- sim.d
      }
      ##### make trees equal with template
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


#' Function for aggregating characters under a specified set of terms (for example, body regions)
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

#' Function to write Ontology CTMC models as a .Rev file for execution in Revbayes
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

