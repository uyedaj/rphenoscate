# FUNCTIONS for to work with Ontology Index

#' @title Convert list to edge matrix
#' @description Takes list of charater annotations amd creates an edge matrix comprising two columns: from and to.
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



######################
#' @title Convert table to list
#' @description Takes table where each row consists of charcter number and ontology annotations and returns a list.
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