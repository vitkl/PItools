##' Filter InteractionSubsetFASTA_list by the size of the interaction sets (number of proteins in each)
##' @name filterInteractionSubsetFASTA_list
##' @author Vitalii Kleshchevnikov
##' @param interactionFASTA_list object of class InteractionSubsetFASTA_list created by \code{\link{listInteractionSubsetFASTA}}
##' @param length_set1_min mininal length of set 1, for example, if a set has less than 1 protein QSLIMFinder wouldn't work.
##' @param length_set2_min mininal length of set 2
##' @param length_set1_max maximal length of set 1, for example, QSLIMFinder has an upper limit of 500 sequences
##' @param length_set2_max maximal length of set 2
##' @return object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @import data.table
##' @import Biostrings
##' @export filterInteractionSubsetFASTA_list
##' @seealso \code{\link{listInteractionSubsetFASTA}}
##' @examples
##' forSLIMFinder_Ready = filterInteractionSubsetFASTA_list(forSLIMFinder,  length_set1_min = 3, length_set2_min = 1)
filterInteractionSubsetFASTA_list = function(interactionFASTA_list, length_set1_min = 0, length_set2_min = 0, length_set1_max = Inf, length_set2_max = Inf) {
  select_both = sapply(1:interactionFASTA_list$length, function(i){
    select_set1 = interactionFASTA_list$interaction_subset[[i]]$length_set1 >= length_set1_min
    select_set2 = interactionFASTA_list$interaction_subset[[i]]$length_set2 >= length_set2_min
    select_set1max = interactionFASTA_list$interaction_subset[[i]]$length_set1 <= length_set1_max
    select_set2max = interactionFASTA_list$interaction_subset[[i]]$length_set2 <= length_set2_max
    select_both = select_set1 & select_set2
    return(select_both)
  })
  interactionFASTA_list$fasta_subset_list = interactionFASTA_list$fasta_subset_list[select_both]
  interactionFASTA_list$interaction_subset = interactionFASTA_list$interaction_subset[select_both]
  interactionFASTA_list$length = length(interactionFASTA_list$fasta_subset_list)
  return(interactionFASTA_list)
}
