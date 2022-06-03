

#' Plots a heatmap along with a phylogeny and trait tree.
#'
#' @export
ontologyHeatMap <- function(td, njt, start=3, margs=c(0.2, 0.25), ...){
  #vals <- na.omit(unique(do.call(c, lapply(3:ncol(td$dat), function(x) unique(as.character(td$dat[[x]]))))))

  X <- do.call(cbind, lapply(start:ncol(td$dat), function(x) as.numeric(recode(td$dat[[x]], "0 and 1"=0.5, "1 and 0"=0.5, "1"=1, "0"=0, "2"=2, "3"=3))))
  colnames(X) <- colnames(td$dat)[start:ncol(td$dat)]

  .vals <- sort(na.omit(unique(as.vector(X))))
  dimx <- dim(X)

  tree1 <- phytools::force.ultrametric(td$phy,method = "extend", message = F)
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

  plot(tree1, x.lim=c(0,(1+margs[2])*(h2$x.lim[2]+h1$x.lim[1])), y.lim=c(0,h1$y.lim[2]+h2$y.lim[2]), show.tip.label = F,...)

  if(!is.null(njt)){
    par(new = TRUE)
    plot(tree2, direction = "downwards", cex = 0.5, x.lim=c(-h1$x.lim[2],h2$x.lim[2]), y.lim=c((-h1$y.lim[2])-0.01*dimx[1],h2$y.lim[2]))
  }

  return ()
}


#' Makes a trait tree using semantic similarity.
#'
#' @export
makeTraitTree <- function (td, external = F, ONT, method = "nj")
{

  if(external == T){

    traits <- colnames(td$dat) # extract trait terms

    # Get subsumer matrix for external ontology
    x <- names(ONT$name[match(traits, ONT$name)])

    anc <- get_ancestors(ONT, x)
    mat <- lapply(as.list(x), function(x) {anc %in% get_ancestors(ONT, x)} )
    mat <- do.call(cbind, mat)

    rownames(mat) <- anc
    colnames(mat) <- x

    smat <- mat*1

    # Get semantic similarity matrix
    semanticSimilarityMatrix <- rphenoscape::tanimoto_similarity(subsumer_mat = smat)
    rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- traits


  }else{

    traits <- colnames(td$dat) # extract trait terms

    # Get semantic similarity matrix
    semanticSimilarityMatrix <- rphenoscape::jaccard_similarity(terms = traits, .colnames = "label", .labels = traits)

    # rownames(semanticSimilarityMatrix) <- colnames(semanticSimilarityMatrix) <- traits

  }

  if(method == "nj"){

    # Neighbor-joining tree of SS matrix
    tree <- nj(1 - semanticSimilarityMatrix)

  }

  if(method == "cls"){

    # Clustering dendrogram of SS matrix
    tree <- hclust(as.dist((1-semanticSimilarityMatrix)))

  }

  return (tree)

}


