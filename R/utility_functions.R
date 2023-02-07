

#####################################
##### UPDATE: Diego (2022)      #####
#####################################

#' Utility function for process an ontoFAST object.
#'
#' @param ontoFAST A named list produced by the function `annot_all_chars` from the package 'ontoFAST'
#' @param ONT An 'ontology_index' object of an external ontology imported in R using the package 'ontologyIndex'
#' @param s.filter Indicates if to use a strict filter of terms
#' @param g.filter Indicates if to use a general filter of terms
#' @param s.terms A character vector with particular ontology terms to be sorted out
#' @param g.terms A character vector with general terms or expressions to be sorted out
#'
#' @return M
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


#' Merge compatible phenotypes.
#'
#' Internal function. Not exported.
#'
merge.pheno <- function(x){

  diag(x) <- 0
  x[x != 1] <- 0

  res <- clusters(graph.adjacency(x))$membership

  return(res)

}


#' Rescore individual states.
#'
#' Internal function. Not exported.
#'
rescore <- function(x,y){

  u <- rep("?", length(x))

  u[x] <- y

  return(u)

}


#' Rescore all states and build characters.
#'
#' Internal function. Not exported.
#'
rescore.all <- function(x, y){

  # Collapse states #
  u <- apply(mapply(x = x, y = y, function(x,y) rescore(x,y)), 1, function(x) paste0(unique(x), collapse = "") )

  # Recoding missings #
  u <- gsub(u, pattern = "\\?", replacement = "")
  u <- gsub(u, pattern = "^$", replacement = "?")

  # Recoding polymorphisms (MrBayes) #
  u <- gsub(u, pattern = "(.)", replacement = "\\1,")
  u <- gsub(u, pattern = ",$", replacement = "")
  u <- gsub(u, pattern = "^(\\d),", replacement = "(\\1,")
  u <- gsub(u, pattern = ",(\\d)$", replacement = ",\\1)")

  return(u)

}


#' Infer clusters of phenotypes based on mutually exclusivity. 
#'
#' @param exclu.obj A list produced by the function `mutually_exclusive` from the package 'rphenoscape'
#'
#' @return A list containing:
#' $phenotypes The inferred clusters of phenotypes (phenotype IDs)
#' $submatrices The exclusivity submatrices corresponding to each phenotype cluster extracted from the $matrix element of the exclu.obj object 
#'
#' \dontrun{
#'     extract.chars(exclu.obj)
#' }
#'
#' @importFrom igraph graph.adjacency
#' @importFrom igraph clusters
#'
#' @export
extract.chars <- function(exclu.obj){

  # Convert exclusivity matrix to binary #
  m <- exclu.obj$matrix

  # Get only mutually exclusive phenotypes #
  #m[m != 5] <- 0
  #m[m == 5] <- 1

  # Get only mutually exclusive and compatible phenotypes #
  m[m != 1 & m != 5] <- 0
  m[m == 1 | m == 5] <- 1

  # Extract all clusters based on mutual exclusivity #
  cls <- clusters(graph.adjacency(m))$membership

  # Extract all exclusivity submatrices #
  S <- split(colnames(exclu.obj$matrix), f = as.factor(cls))
  S <- lapply(S, function(x) exclu.obj$matrix[c(x),c(x)] )

  # Define clusters of phenotypes #
  x <- unique(c(exclu.obj$dataframe$id.1, exclu.obj$dataframe$id.2))
  CH <- split(x, f = as.factor(cls))

  # Return results #
  return(list(phenotypes = CH, submatrices = S))

}


#' Build characters based on inferred phenotype clusters.
#'
#' @param CH A named list produced by the function `extract.chars`
#' @param chars.obj A data.frame or list produced by the function `chars` from the package 'rphenoscape'
#'
#' @return A list with two main elements: $solved and $unsolved
#' $solved contains all phenotypes assigned to inferred clusters
#' $unsolved contains phenotypes that were not assigned to clusters or singletons. Each main element contains:
#' $chars Inferred clusters of phenotypes with character statements writen in natural language, as in the original study
#' $clusters Inferred clusters shown as phenotype IDs
#' $tokens Codes assigned to character states
#'
#' \dontrun{
#'     build.chars(CH, chars.obj)
#' }
#'
#' @export
build.chars <- function(CH, chars.obj){

  ## STEP 1: Isolate different types of phenotype clusters ##
  # Extract 'singletons' #
  CH.sing <- CH$phenotypes[sapply(CH$phenotypes, function(x) length(x) < 2 )]

  # Extract clusters with multiple phenotypes #
  CH.mult <- CH$phenotypes[!sapply(CH$phenotypes, function(x) length(x) < 2 )]

  # Extract submatrices for clusters with multiple phenotypes #
  S <- CH$submatrices[names(CH.mult)]

  # Extract submatrices with mutually compatible phenotypes (= should merge) #
  S <- S[sapply(S, function(x) any(x[upper.tri(x)] != 5) )]

  # Extract clusters with mutually exclusive and compatible phenotypes #
  CH.mult.merge <- CH.mult[names(S)]

  # Extract clusters with only mutually exclusive phenotypes #
  CH.mult.done <- CH.mult[!names(CH.mult) %in% names(S)]



  ## STEP 2: Deal with mutually compatible phenotypes ##
  # Get phenotypes that should be merged #
  m <- lapply(S, merge.pheno)

  # Get initial clusters of phenotypes to be merged #
  M <- mapply(x = CH.mult.merge, y = m, function(x,y) split(x, f = as.factor(y)), SIMPLIFY = FALSE )

  # Match phenotype clusters with state labels in the original works to check for redundant information #
  x <- sapply(M, function(x) any(sapply(x, function(x) length(unique(chars.obj$state.label[match(x, chars.obj$phenotype.id)])) > 1 )) )

  # Extract clusters of phenotype that can be merged #
  M.merge <- M[!x]

  # Extract clusters of phenotype that cannot be merged immediately #
  M.redo <- M[x]

  ## STEP 3: Get character and character state labels ##

  ## Type 1: no problems ##
  # Copy original phenotypes #
  chars1 <- CH.mult.done

  # Get character labels for 'regular' phenotypes (only mutually exclusive) #
  char.lab1 <- lapply(chars1, function(x) unique(chars.obj$character.label[match(x, chars.obj$phenotype.id)]) )
  names(chars1) <- sapply(char.lab1, function(x) paste0(x, collapse = " + ") )

  # Get state labels for 'regular' phenotypes (only mutually exclusive) #
  chars1 <- lapply(chars1, function(x) chars.obj$state.label[match(x, chars.obj$phenotype.id)] )

  ## Type 2: compatible phenotypes merged and matching original state labels ##
  # Copy original phenotypes #
  chars2 <- M.merge

  # Get character labels for mutually compatible phenotypes #
  char.lab2 <- lapply(chars2, function(x) sapply(x, function(x)
    unique(chars.obj$character.label[match(x, chars.obj$phenotype.id)]) ) )
  char.lab2 <- lapply(char.lab2, function(x) paste0(unique(x), collapse = " + ") )
  names(chars2) <- char.lab2

  # Get state labels for mutually compatible phenotypes #
  chars2 <- lapply(chars2, function(x) sapply(x, function(x)
    paste0(unique(chars.obj$state.label[match(x, chars.obj$phenotype.id)]),collapse = " + ") ) )
  chars2 <- lapply(chars2, unname)

  if(length(M.redo) > 0){

    ## STEP 2: Deal with mutually compatible phenotypes ##

    # Extract cluster of phenotypes that require merge of non-redundant state labels #
    M.redo.mult <- M.redo[sapply(m[names(M.redo)], function(x) length(unique(x)) > 1 )]

    # Extract clusters with 'singletons' (in this case, all phenotypes merging to only one 'state') #
    M.redo.sing <- M.redo[sapply(m[names(M.redo)], function(x) length(unique(x)) == 1 )]

    # Join unsolved phenotypes #
    M.redo.unsolved <- c(CH.sing, M.redo.sing)

    ## STEP 3: Get character and character state labels ##

    ## Type 3: compatible phenotypes merged and non-matching original state labels (= semantic diversity) ##
    # Copy original phenotypes #
    chars3 <- M.redo.mult

    # Get character labels for mutually compatible phenotypes with non-redundant mergeable state labels #
    char.lab3 <- lapply(chars3, function(x) sapply(x, function(x)
      unique(chars.obj$character.label[match(x, chars.obj$phenotype.id)]) ) )
    char.lab3 <- lapply(char.lab3, function(x) paste0(unique(unlist(x)), collapse = " + ") )
    names(chars3) <- char.lab3

    # Get state labels for mutually compatible phenotypes with non-redundant mergeable state labels #
    chars3 <- lapply(chars3, function(x)
      unname(sapply(x, function(x) unname(paste0(unique(chars.obj$state.label[match(x, chars.obj$phenotype.id)]),collapse = " + ")) )) )

    ## Type 4: non-solved: phenotypes not clustering or all clustering as a single compatible 'state'  ##
    # Copy original phenotypes #
    chars4 <- M.redo.unsolved

    # Get character labels for new clusters of phenotypes from strategy 2 #
    char.lab4 <- lapply(chars4, function(x) sapply(x, function(x)
      unique(chars.obj$character.label[match(x, chars.obj$phenotype.id)]) ) )
    char.lab4 <- lapply(char.lab4, function(x) paste0(unique(unlist(x)), collapse = " + ") )
    names(chars4) <- char.lab4

    # Get state labels for new clusters of phenotypes from strategy 2 #
    chars4 <- lapply(chars4, function(x)
      unname(lapply(x, function(x) unname(unique(chars.obj$state.label[match(x, chars.obj$phenotype.id)])) )) )
    chars4 <- lapply(chars4, unlist)

    ## STEP 4: Build characters ##
    # Create a list to store results #
    res <- list()

    # Join character statements #
    res$solved$chars <- c(chars1, chars2, chars3)
    res$unsolved$chars <- chars4

    # Join clusters of phenotypes $
    res$solved$clusters <- c(CH.mult.done, M.merge, M.redo.mult)
    res$unsolved$clusters <- M.redo.unsolved

    # Set character state tokens #
    res$solved$tokens <- lapply(res$solved$chars, function(x) paste(0:(length(x) - 1)) )
    res$unsolved$tokens <- lapply(res$unsolved$chars, function(x) paste(0:(length(x) - 1)) )

  }else{

    ## STEP 4: Build characters ##
    # Create a list to store results #
    res <- list()

    # Join character statements #
    res$solved$chars <- c(chars1, chars2)

    # Join clusters of phenotypes $
    res$solved$clusters <- c(CH.mult.done, M.merge)

    # Set character state tokens #
    res$solved$tokens <- lapply(res$solved$chars, function(x) paste(0:(length(x) - 1)) )

  }

  # Return results #
  return(res)

}


#' Build character matrices based on inferred characters.
#'
#' @param character A character vector with the taxon names to be scored
#' @param pheno.obj A list of phenotypes produced by the function `get_phenotypes` and converted with `as.phenotype` from the package 'rphenoscape'
#' @param chars.obj A data.frame or list produced by the function `chars` from the package 'rphenoscape'
#' @param CHARS A named list of inferred characters produced by the function `build.chars`
#'
#' \dontrun{
#'     build.matrix(tax, pheno.obj, chars.obj, CHARS)
#' }
#'
#' @export
build.matrix <- function(tax, pheno.obj, chars.obj, CHARS){

  # Get set of taxa for each phenotype #
  tax.list <- lapply(CHARS$solved$clusters, function(y)
    lapply(y, function(x) unlist(sapply(pheno.obj, function(x)
      x$taxa$label)[match(unlist(x), sapply(pheno.obj, function(x) x$id ))]) ))

  # Check for taxa exhibiting each phenotype #
  tax.scores <- lapply(tax.list, function(x) lapply(x, function(x) tax %in% x) )

  # Get character matrix #
  M <- mapply(x = tax.scores, y = CHARS$solved$tokens, function(x,y) rescore.all(x, y), SIMPLIFY = FALSE )

  # Build character matrix #
  M <- as.data.frame(do.call(cbind, M))
  rownames(M) <- tax
  colnames(M) <- paste0("C", 1:length(tax.list))

  # Return results #
  return(M)

}


#########################
##### OLD FUNCTIONS #####
#########################

#' Utility function for displaying the number of matching species and average value across them from an Ontotrace matrix.
#'
#' @param x Ontotrace matrix
#'
#' @return Prints the number of matching species and average value across them from an Ontotrace
#'
#' @importFrom stats na.omit
#'
#' @export
print_coverage <- function(x){
  coverage <- apply(x, 2, function(x) sum(!is.na(x)))
  percentage <- round(coverage/dim(x)[1],2)
  average <- sapply(x, function(x) round(mean(as.numeric(na.omit(x[x %in% c("0", "1")])), na.rm = TRUE),2))
  cover <- dplyr::tibble(traits = names(coverage), coverage, percentage, average)
  tmp <- dplyr::filter(cover, coverage > 0, average < 1, average > 0) %>% arrange(desc(coverage))
  print(tmp)
}


#' Utility function for filtering data based on missing traits and taxa.
#'
#' @param td A treedata object to filter
#' @param traits The traits to filter
#' @param taxa The taxa to filter
#'
#' @return td The filtered treedata object
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
#' @param char_info The character data table to clean
#' @param dep.mat Dependency matrix
#' @param td The treedata object with amalgamated characters
#'
#' @return char_info_comb Cleaned character data table after amalgamating characters.
#'
#' @importFrom stats setNames
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
#' @param x A full IRI
#'
#' @return x The IRI with its url stripped.
#'
#' @export
strip_IRI <- function(x){
  x <- gsub("http://purl.obolibrary.org/obo/", "", x)
  x <- gsub("_", ":", x)
  return(x)
}
