##' Recode names of the AAStringSet object
##' @name recodeFASTA
##' @author Vitalii Kleshchevnikov
##' @param fasta AAStringSet object for which to generate new names or the interactionSubsetFASTA which to rename back to the original names
##' @param mapping data.table, mapping table between new and original names
##' @param i which element of the interactionSubsetFASTA to rename back to the original names
##' @return list of the original AAStringSet object with new names and mapping table between new and original names
##' @return interactionSubsetFASTA object with original names as provided by mapping table
##' @import Biostrings
##' @import data.table
##' @export recodeFASTA
##' @seealso \code{\link{downloadFasta_postproc}}
##' @examples
##' # from original to new
##' recoded.fasta = recodeFASTA(fasta)
##' # from new to original
##' subset_fasta = recodeFASTA(subset_fasta, recoded.fasta$names_mapping, i = 1)
recodeFASTA = function(fasta, mapping = NULL, i = NULL){
  if(is.null(mapping)){
    old_names = names(fasta)
    new_names = paste0("P",seq_along(old_names))
    names(fasta) = new_names
    names_mapping = data.table(old_names = old_names, new_names = new_names)
    return(list(fasta = fasta, names_mapping = names_mapping))
  } else {
    if(class(fasta) == "interactionSubsetFASTA"){
      if(!is.null(i)){
        names(fasta$fasta_subset_list[[i]]) = mapping$old_names[match(names(fasta$fasta_subset_list[[i]]), mapping$new_names)]
        return(fasta)
      } else stop("i not provided")
    } else stop("fasta is not interactionSubsetFASTA")
  }
}
