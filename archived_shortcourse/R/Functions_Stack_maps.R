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


#### stack two discrete stm's lists; x,y are the list of state names (i.e. maps)
stack2<-function(x,y){
  mapply(function(x,y) 
  {paste(x,y, sep="") },
  x=x, y=y )
}

#' Final stack of maps for a set of stochastic maps stored in a directory
#' cc chars id to stack
#' ntrees number of trees to stack
#' dirW directory for zip file
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


#' Final stack of maps for a set of stochastic maps stored in a list
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
