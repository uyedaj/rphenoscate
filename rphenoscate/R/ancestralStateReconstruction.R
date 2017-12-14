library(treeplyr)
## Read in data files
td <- readRDS("../data/matchedTreeDataDated.rds")
td <- select(td, -1)
## Read in relation table
rt <- read.csv("../rphenoscate/R/dependencies_with_names.txt")
rt <- apply(rt, 2, function(x) gsub(" ", ".", x, fixed=TRUE))
rt <- apply(rt, 2, function(x) gsub("-", ".", x, fixed=TRUE))

## Find variable dependent traits for analysis
variableDep <- NULL
tds <- list()
for(i in 1:nrow(rt)){
  t1 <- rt[i,2]
  t2 <- rt[i,4]
  tmptd <- select(td, t1, t2) %>% filter(., (!is.na(td[[t1]])), (!is.na(td[[t2]]))) %>% mutate(., X = paste0(.[[2]], .[[1]]))
  tmpvar <- apply(tmptd$dat, 2, function(x) length(unique(x)))
  if(all(tmpvar > 1) & tmpvar[3] > 2) variableDep <- c(variableDep, 1) else { variableDep <- c(variableDep, 0) }
  tds <- c(tds, list(tmptd))
}

dependentCases <- which(variableDep==1)
for(i in dependentCases) print(table(tds[[i]][['X']]))

pdf("../data/asrDependentCases.pdf")
for(i in dependentCases)
{
    tmptd <- tds[[i]]
    phy <- multi2di(tmptd$phy, random=FALSE)
    phy$edge.length[phy$edge.length==0] <- .Machine$double.eps
    aceDep <- ace(tmptd[['X']], phy , "discrete", model="ARD", marginal=TRUE)
    aceInd1 <- ace(tmptd[[1]], phy, "discrete", model="ARD", marginal=TRUE)
    aceInd2 <- ace(tmptd[[2]], phy, "discrete", model="ARD", marginal=TRUE)

    adjxy <- c(3,3)
    par(mfrow=c(1,2))
    plot(phy, show.tip.label=FALSE)
                                        #nodelabels(pie=aceDep$lik.anc, piecol=c("black", "green", "red"))
    try(nodelabels(pie=cbind(aceDep$lik.anc[,1], aceDep$lik.anc[,2]+aceDep$lik.anc[,3]), piecol=c("black", "green"), adj=-1*adjxy))
    try(nodelabels(pie=cbind(aceDep$lik.anc[,1] + aceDep$lik.anc[,2], aceDep$lik.anc[,3]), piecol=c("black", "red"), adj=adjxy))
    tiplabels(pch=21, col=tmptd[[1]]+1, bg=tmptd[[1]]+1, cex=0.5)
    tiplabels(pch=21, bg=c("black", "green")[tmptd[[2]]+1], col=c("black", "green")[tmptd[[2]]+1],, adj=c(7.5,0.5), cex=0.5)

    plot(phy, show.tip.label=FALSE)
    try(nodelabels(pie=aceInd1$lik.anc, piecol=c("black", "red"), adj=adjxy))
    try(nodelabels(pie=aceInd2$lik.anc, piecol=c("black", "green"), adj=-1*adjxy))
    tiplabels(pch=21, col=tmptd[[1]]+1, bg=tmptd[[1]]+1, cex=0.5)
    tiplabels(pch=21, bg=c("black", "green")[tmptd[[2]]+1], col=c("black", "green")[tmptd[[2]]+1],, adj=c(7.5,0.5), cex=0.5)
}
dev.off()

## Phytools
library(phytools)
aceDep.simmap <- make.simmap(phy, tmptd[['X']],nsim=10, model="ARD")
aceInd.simmap <- make.simmap(phy, tmptd[[1]],nsim=10, model="ARD")

plotSimmap(aceDep.simmap[[1]], colors = setNames(c("black", "green", "red"), c("00", "10", "11")), ftype="off")
tiplabels(pch=21, bg=tmptd[[1]]+1)
tiplabels(pch=21, bg=c("black", "green")[tmptd[[2]]+1], adj=c(5,0.5))
plotSimmap(aceInd.simmap[[1]], ftype="off")
tiplabels(pch=21, bg=tmptd[[1]]+1)
tiplabels(pch=21, bg=c("black", "green")[tmptd[[2]]+1], adj=c(5,0.5))

