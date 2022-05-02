
########################################
##### Stochastic Mapping functions #####
########################################

## Imported from OntoPhylo ##

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


## Imported from OntoPhylo ##

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


## Imported from OntoPhylo ##

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


## Imported from OntoPhylo ##

#' Function for processing corHMM stochastic maps.
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

## Imported from OntoPhylo ##

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

## Imported from OntoPhylo ##

## Not exported ##
# Stack two discrete stm's lists; x,y are the list of state names (i.e. maps).
stack2<-function(x,y){
  mapply(function(x,y)
  {paste(x,y, sep="") },
  x=x, y=y )
}
