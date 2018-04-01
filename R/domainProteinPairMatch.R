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

  # for testing
  #envir = new.env()
  #load("../viral_project/processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Full_IntAct3.FALSE_clust201802.RData", envir = envir)
  #domain_res_env = R.utils::env(load("../viral_project/processed_data_files/what_we_find_VS_ELM_clust20171019.RData"))
  #domain_res = domain_res_env$res_count
  #dbfile_main = "../viral_project/data_files/instances_all.gff"
  #dburl_main = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic="
  #instances_all = ELMdb2GRanges(dbfile = dbfile_main,
  #                              dburl = dburl_main,
  #                              tsvurl = gsub("gff", "tsv", dburl_main),
  #                              tsvfile = gsub("gff", "tsv", dbfile_main))
  #grange = instances_all
  #interactionSubsetFASTA = envir$forSLIMFinder_Ready
  # end for testing

  interaction_subset = InteractionSubsetFASTA_list$interaction_subset
  nodeX = formula(paste0("~",domain_res$nodes$nodeX))[[2]]
  nodeY = formula(paste0("~",domain_res$nodes$nodeY))[[2]]
  which_to_keep = sapply(interaction_subset, function(interaction_subset_x, domain_res, nodeX, nodeY){
    name = interaction_subset_x$name
    name = unlist(strsplit(name, ":"))
    which_nodeY = domain_res$data_with_pval[, eval(nodeY)] %in% name[1]
    return(domain_res$data_with_pval[which_nodeY, name[2] %in% eval(nodeX)])
  }, domain_res, nodeX, nodeY)
  if(remove){
    InteractionSubsetFASTA_list$fasta_subset_list = InteractionSubsetFASTA_list$fasta_subset_list[which_to_keep]
    InteractionSubsetFASTA_list$interaction_subset = InteractionSubsetFASTA_list$interaction_subset[which_to_keep]
    InteractionSubsetFASTA_list$length = length(InteractionSubsetFASTA_list$fasta_subset_list)
    return(InteractionSubsetFASTA_list)
  } else return(which_to_keep)
}
