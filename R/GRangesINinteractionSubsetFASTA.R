##' Check if proteins in GRanges were used for motif search using QSLIMFinder (SLIMFinder)
##' @rdname GRangesINinteractionSubsetFASTA
##' @name GRangesINinteractionSubsetFASTA
##' @author Vitalii Kleshchevnikov
##' @param GRanges object containing instances of linear motifs or domains in proteins (identified by UniProtKB accession)
##' @param file object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @return list containing 2 vectors: which GRanges are in InteractionSubsetFASTA_list and seqlengths vector for these GRanges
##' @import data.table
##' @import GenomicRanges
##' @export GRangesINinteractionSubsetFASTA
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{runQSLIMFinder}}, \code{\link{SLIMFinderOcc2GRanges}}, \code{\link{listInteractionSubsetFASTA}}, \code{\link{randomGRanges}}
GRangesINinteractionSubsetFASTA = function(grange, interactionSubsetFASTA) {
  sequencesSearched = sapply(1:length(interactionSubsetFASTA$fasta_subset_list), function(i) {
    x = width(interactionSubsetFASTA$fasta_subset_list[[i]])
    names(x) = names(interactionSubsetFASTA$fasta_subset_list[[i]])
    x
    })
  sequencesSearched = Reduce(c, sequencesSearched)
  sequencesSearched = sequencesSearched[unique(names(sequencesSearched))]

  list(granges_in_sequencesSearched = as.character(seqnames(grange)) %in% names(sequencesSearched), seqlengths = sequencesSearched)
}
