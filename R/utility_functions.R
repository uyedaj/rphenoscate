

#' Utility function for displaying the number of matching species and average value across them from an Ontotrace matrix.
#'
#' @export
print_coverage <- function(x){
  coverage <- apply(x, 2, function(x) sum(!is.na(x)))
  average <- sapply(x, function(x) mean(as.numeric(na.omit(x[x %in% c("0", "1")])), na.rm=TRUE))
  cover <- cbind(coverage, average)
  tmp <- dplyr::filter(data.frame(traits=rownames(cover), cover), coverage > 0, average < 1, average > 0) %>% arrange(., desc(coverage))
  print(tmp)
}


#' Utility function for filtering data based on missing traits and taxa.
#'
#' @export
filter_coverage <- function(td, traits=0, taxa=0){
  taxa_coverage <- apply(td$dat, 1, function(x) mean(as.numeric(!is.na(x))))
  trait_coverage <- apply(td$dat, 2, function(x) mean(as.numeric(!is.na(x))))
  td <- dplyr::filter(td, taxa_coverage > taxa)
  td <- dplyr::select(td, which(trait_coverage > traits))
  return(td)
}


#' Utility function for cleaning up character data table after amalgamating characters.
#'
#' @export
dropDependentTraits <- function(char_info, dep.mat, td){
  char_info_comb <- char_info[which(apply(dep.mat, 1, sum, na.rm=TRUE)==0), c(1,5)]
  new.traits <- colnames(td$dat)
  old.traits <- sapply(new.traits, function(x) strsplit(x, "+", fixed=TRUE)[[1]][1])
  trait.trans <- setNames(new.traits, old.traits)
  char_info_comb$ID <- unname(trait.trans[as.character(char_info_comb$ID)])
  return(char_info_comb)
}


#' Utility function for stripping off url from IRI.
#'
#' @export
strip_IRI <- function(x){
  x <- gsub("http://purl.obolibrary.org/obo/", "", x)
  x <- gsub("_", ":", x)
  return(x)
}
