##' Check if proteins in GRanges were used for motif search using QSLIMFinder (SLIMFinder)
##' @rdname GRangesINinteractionSubsetFASTA
##' @name GRangesINinteractionSubsetFASTA
##' @author Vitalii Kleshchevnikov
##' @param GRanges object containing instances of linear motifs or domains in proteins (identified by UniProtKB accession)
##' @param interactionSubsetFASTA object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @param query_only return only GRanges for query proteins
##' @param clustermq Use clustermq LSF job scheduler (TRUE) as an alternative to parLapply (FALSE). Details: \link[clustermq]{Q}
##' @param clustermq_seed When using clustermq: Seed for random number generation.
##' @param clustermq_memory When using clustermq: memory requested for each job
##' @param clustermq_job_size When using clustermq: The number of function calls per job
##' @return list containing 2 vectors: which GRanges are in InteractionSubsetFASTA_list and seqlengths vector for these GRanges
##' @import data.table
##' @import GenomicRanges
##' @export GRangesINinteractionSubsetFASTA
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{runQSLIMFinder}}, \code{\link{SLIMFinderOcc2GRanges}}, \code{\link{listInteractionSubsetFASTA}}, \code{\link{randomGRanges}}
GRangesINinteractionSubsetFASTA = function(grange, interactionSubsetFASTA, query_only = F,
                                           clustermq = F, clustermq_seed = 128965,
                                           clustermq_memory = 2000, clustermq_job_size = 100) {
  if(class(interactionSubsetFASTA) != "InteractionSubsetFASTA_list") stop("interactionSubsetFASTA should be of class InteractionSubsetFASTA_list returned by listInteractionSubsetFASTA()")
  fasta_subset_list = interactionSubsetFASTA$fasta_subset_list
  sequences_Searched = sapply(fasta_subset_list, function(fasta_subset) {
    x = width(fasta_subset)
    names(x) = names(fasta_subset)
    x
  })
  sequences_Searched = Reduce(c, sequences_Searched)
  sequences_Searched = sequences_Searched[unique(names(sequences_Searched))]

  if(query_only){
    interaction_subset = interactionSubsetFASTA$interaction_subset
    query_Names = sapply(interaction_subset, function(interaction_sub) {
      interaction_sub$ids_set2
    })

    sequences_Searched = sequences_Searched[unique(query_Names)]
    sequences_Searched = sequences_Searched[unique(names(sequences_Searched))]
  }

  list(granges_in_sequences_Searched = as.character(seqnames(grange)) %in% names(sequences_Searched), seqlengths = sequences_Searched)
}
