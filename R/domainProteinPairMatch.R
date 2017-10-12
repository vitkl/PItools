##' Check if the protein in the set2 is associated with domain that is contained in the protein in the name and remove if not
##' @name domainProteinPairMatch
##' @author Vitalii Kleshchevnikov
##' @param InteractionSubsetFASTA_list object of class InteractionSubsetFASTA_list created by \code{\link{listInteractionSubsetFASTA}}
##' @param domain_res the result of domain enrichment analysis: object class XYZinteration_XZEmpiricalPval, the output of \code{\link{permutationPval}}. The domain should be in nodeZ, the set2 protein should be in the nodeX, the set1 protein should be in the nodeY.
##' @param remove logical, remove and return InteractionSubsetFASTA_list or just show which sequence and interaction sets have a protein in the set2 associated with domain that is contained in the protein in the name.
##' @return object of class InteractionSubsetFASTA_list filtered for matching protein and containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @import data.table
##' @import Biostrings
##' @export domainProteinPairMatch
##' @seealso \code{\link{listInteractionSubsetFASTA}}, \code{\link{permutationPval}}
##' @examples
##' forSLIMFinder_Ready = domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = T)
domainProteinPairMatch = function(InteractionSubsetFASTA_list, domain_res, remove = T) {
  which_to_keep = sapply(1:InteractionSubsetFASTA_list$length, function(i, InteractionSubsetFASTA_list, domain_res){
    name = InteractionSubsetFASTA_list$interaction_subset[[i]]$name
    name = unlist(strsplit(name, ":"))
    nodeX = formula(paste0("~",domain_res$nodes$nodeX))[[2]]
    nodeY = formula(paste0("~",domain_res$nodes$nodeY))[[2]]
    nodeZ = formula(paste0("~",domain_res$nodes$nodeZ))[[2]]
    which_nodeY = domain_res$data_with_pval[, eval(nodeY)] %in% name[1]
    return(domain_res$data_with_pval[which_nodeY, name[2] %in% eval(nodeX)])
  }, InteractionSubsetFASTA_list, domain_res)

  if(remove){
    InteractionSubsetFASTA_list$fasta_subset_list = InteractionSubsetFASTA_list$fasta_subset_list[which_to_keep]
    InteractionSubsetFASTA_list$interaction_subset = InteractionSubsetFASTA_list$interaction_subset[which_to_keep]
    InteractionSubsetFASTA_list$length = length(InteractionSubsetFASTA_list$fasta_subset_list)
    return(InteractionSubsetFASTA_list)
  } else return(which_to_keep)
}
