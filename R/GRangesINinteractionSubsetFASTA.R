##' Check if proteins in GRanges were used for motif search using QSLIMFinder (SLIMFinder)
##' @rdname GRangesINinteractionSubsetFASTA
##' @name GRangesINinteractionSubsetFASTA
##' @author Vitalii Kleshchevnikov
##' @param GRanges object containing instances of linear motifs or domains in proteins (identified by UniProtKB accession)
##' @param interactionSubsetFASTA object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @return list containing 2 vectors: which GRanges are in InteractionSubsetFASTA_list and seqlengths vector for these GRanges
##' @import data.table
##' @import GenomicRanges
##' @export GRangesINinteractionSubsetFASTA
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{runQSLIMFinder}}, \code{\link{SLIMFinderOcc2GRanges}}, \code{\link{listInteractionSubsetFASTA}}, \code{\link{randomGRanges}}
GRangesINinteractionSubsetFASTA = function(grange, interactionSubsetFASTA, query_only = F) {
  if(class(interactionSubsetFASTA) != "InteractionSubsetFASTA_list") stop("interactionSubsetFASTA should be of class InteractionSubsetFASTA_list returned by listInteractionSubsetFASTA()")

  sequences_Searched = sapply(1:length(interactionSubsetFASTA$fasta_subset_list), function(i) {
    x = width(interactionSubsetFASTA$fasta_subset_list[[i]])
    names(x) = names(interactionSubsetFASTA$fasta_subset_list[[i]])
    x
    })
  sequences_Searched = Reduce(c, sequences_Searched)
  sequences_Searched = sequences_Searched[unique(names(sequences_Searched))]

  if(query_only){
    query_Names = sapply(1:length(interactionSubsetFASTA$interaction_subset), function(i) {
      interactionSubsetFASTA$interaction_subset[[i]]$ids_set2
    })
    sequences_Searched = sequences_Searched[query_Names]
    sequences_Searched = sequences_Searched[unique(names(sequences_Searched))]
  }

  list(granges_in_sequences_Searched = as.character(seqnames(grange)) %in% names(sequences_Searched), seqlengths = sequences_Searched)
}
