#####################################
##### UPDATE: Diego (2022)      #####
#####################################

#' Utility function for process an ontoFAST object
#'
#' @param ontoFAST A named list produced by the function `annot_all_chars` from the package ontoFAST
#' @param ONTO An 'ontology_index' object of an external ontology imported in R using the package 'ontologyIndex'
#' @param s.filter Indicates if to use a strict filter of terms
#' @param g.filter Indicates if to use a general filter of terms
#' @param s.terms A character vector with particular ontology terms to be sorted out
#' @param g.terms A character vector with general terms or expressions to be sorted out
#'
#' @export
process.ontofast <- function(ontoFAST, ONT, s.filter = F, g.filter = F, s.terms, g.terms){

  M <- lapply(ontoFAST, function(x) ONT$name[names(ONT$name) %in% x] )

  if(s.filter == T){

    M <- lapply(M, function(x) x[!x %in% s.terms])

  }

  if(g.filter == T){

    for(j in 1:length(g.terms)){

      M <- lapply(M, function(x) x[!grepl(x, pattern = g.terms[j])] )

    }

  }

  for(i in 1:length(M)){

    print(paste0("CHAR_", i, ": ", names(M)[i]))
    print(rbind(1:length(M[[i]]), M[[i]]))

    x <- as.numeric(readline("Choose a number or 0 to skip: "))

    if(x > 0){ M[[i]] <- M[[i]][x] }else{

      y <- readline("Do you want to try to find a more suitable term? yes or no: ")

      if(y == "yes"){

        z <- as.character(readline("Please, choose a term from your selected ontology (no quotes!): "))

        candidates <- ONT$name[grepl(ONT$name, pattern = z)]

        if(length(candidates) > 1){

          print(cbind(candidates, 1:length(candidates)))

          w <- as.numeric(readline("Please, select one of the terms from the list: "))

          M[[i]] <- candidates[w]

          names(M)[i] <- candidates[w]

        }else{

          M[[i]] <- candidates
          names(M)[i] <- candidates

          }


      }

    }

  }

  return(M)

}


#########################
##### OLD FUNCTIONS #####
#########################

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
