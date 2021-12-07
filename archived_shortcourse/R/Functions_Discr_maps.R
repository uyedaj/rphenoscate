#' make branch length and maps equal, use tree.tmp with rounded (3) edge.length
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
########


#' Reading unsummarized Stoch Map files from ReVBayes
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
#########


#' Reading unsummarized Stoch Map files from ReVBayes for one tree
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
    ######
    #sprintf("%.54f", c(TR[length(TR)], YY[nrow(YY), 2]) )
    # TR[1]==TR[2]
    # all.equal(TR[1], TR[2])
    # all.equal(TR)
    # duplicated(TR)
    #length(int.out)
    #all.equal(TR[length(TR)], YY[nrow(YY), 2])
    #TR[length(TR)]==YY[nrow(YY), 2]
    ######
    
    #TR[length(TR)]==YY[nrow(YY), 2]
    int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE, all.inside = TRUE)
    #int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE)
    #findInterval(seq(0.1, 4, .1), c(0, 0.5, 0.7, 1.5, 1.6, 4 ), left.open=T, rightmost.closed = F)
    maps.n[[i]]<-setNames(YY[,2]-YY[,1], names(tree$maps[[i]])[int.out])
  }
  tree$maps<-maps.n
  return(tree)        
}


#' Reading unsummarized Stoch Map files from ReVBayes
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



##################################################
