##' Check if the protein in the set2 is associated with domain that is contained in the protein in the name and remove if not
##' @name domainProteinPairMatch
##' @author Vitalii Kleshchevnikov
##' @param InteractionSubsetFASTA_list object of class InteractionSubsetFASTA_list created by \code{\link{listInteractionSubsetFASTA}}
##' @param domain_res the result of domain enrichment analysis: object class XYZinteration_XZEmpiricalPval, the output of \code{\link{permutationPval}}. The domain should be in nodeZ, the set2 protein should be in the nodeX, the set1 protein should be in the nodeY.
##' @param non_query_domain_res the result of domain enrichment analysis for non-query proteins, when provided will be used for filtering.
##' @param non_query_domains_N the number of non-query proteins with predicted domains for each dataset. Used only when non_query_domain_res is not NULL
##' @param non_query_set_only If TRUE interacting partners of a seed are taken only from non_query_domain_res, if FALSE - both from non_query_domain_res and domain_res. Used only when non_query_domain_res is not NULL
##' @param query_domains_only If TRUE interacting partners of a seed that are taken from non_query_domain_res must be predicted to interact with a seed using domains predicted for query. Used only when non_query_domain_res is not NULL
##' @param remove logical, remove and return InteractionSubsetFASTA_list or just show which sequence and interaction sets have a protein in the set2 associated with domain that is contained in the protein in the name.
##' @return object of class InteractionSubsetFASTA_list filtered for matching protein and containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @import data.table
##' @import Biostrings
##' @export domainProteinPairMatch
##' @seealso \code{\link{listInteractionSubsetFASTA}}, \code{\link{permutationPval}}
##' @examples
##' forSLIMFinder_Ready = domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = T)
domainProteinPairMatch = function(InteractionSubsetFASTA_list, domain_res, non_query_domain_res = NULL, non_query_domains_N = 0, non_query_set_only = T, query_domains_only = F, remove = T) {

  # for testing
  #envir = R.utils::env(load("../viral_project/processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Full_IntAct3.FALSE_clust201802.RData"))
  #domain_res_env = R.utils::env(load("../viral_project/processed_data_files/what_we_find_VS_ELM_clust20171019.RData"))
  #domain_res = domain_res_env$res_count
  #non_query_domain_res_env = R.utils::env(load("../viral_project/processed_data_files/predict_domain_human_clust20180817.RData"))
  #non_query_domain_res = non_query_domain_res_env$res_count_all
  #non_query_domain_res$data_with_pval = non_query_domain_res$data_with_pval[p.value < 0.5]
  #dbfile_main = "../viral_project/data_files/instances_all.gff"
  #dburl_main = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic="
  #instances_all = ELMdb2GRanges(dbfile = dbfile_main,
  #                              dburl = dburl_main,
  #                              tsvurl = gsub("gff", "tsv", dburl_main),
  #                              tsvfile = gsub("gff", "tsv", dbfile_main))
  #grange = instances_all
  #InteractionSubsetFASTA_list = envir$forSLIMFinder_Ready
  # end for testing

  interaction_subset = InteractionSubsetFASTA_list$interaction_subset
  nodeX = formula(paste0("~",domain_res$nodes$nodeX))[[2]]
  nodeY = formula(paste0("~",domain_res$nodes$nodeY))[[2]]
  nodeZ = formula(paste0("~",domain_res$nodes$nodeZ))[[2]]
  if(!is.null(non_query_domain_res)){
    nodeX_non_query_domain_res = formula(paste0("~",non_query_domain_res$nodes$nodeX))[[2]]
    nodeY_non_query_domain_res = formula(paste0("~",non_query_domain_res$nodes$nodeY))[[2]]
    nodeZ_non_query_domain_res = formula(paste0("~",non_query_domain_res$nodes$nodeZ))[[2]]
  }

  # create variable to be filled
  which_to_keep = logical(length(interaction_subset))
  seq2keep = data.table(seq2keep = character(), name = character())
  # set up progress bar
  pb_check <- progress::progress_bar$new(
    format = "filtering datasets by enriched domains [:bar] :current/:total eta: :eta",
    total = length(interaction_subset), clear = FALSE, width= 80, show_after = 0)
  for (i in 1:length(interaction_subset)) {
    pb_check$tick()
    interaction_subset_x = interaction_subset[[i]]
    name = interaction_subset_x$name
    name = unlist(strsplit(name, ":"))
    which_nodeY = domain_res$data_with_pval[, eval(nodeY)] == name[1]
    if(sum(which_nodeY) >= 1){
      to_keep = domain_res$data_with_pval[which_nodeY, name[2] == eval(nodeX)]
      to_keep = isTRUE(sum(to_keep) >= 1 & is.logical(to_keep))
    } else to_keep = FALSE


    # filter non-query proteins by domain
    if(!is.null(non_query_domain_res)){
      all_ids_set1 = interaction_subset_x$ids_set1
      non_query_which_nodeY = non_query_domain_res$data_with_pval[, eval(nodeY_non_query_domain_res)] == name[1]
      # if non-query domains should be the same as query domains - do filter by this criteria
      if(query_domains_only){
        query_domains = domain_res$data_with_pval[which_nodeY, eval(nodeZ)]
        which_set1_proteins = non_query_domain_res$data_with_pval[non_query_which_nodeY &
                                                                    eval(nodeZ_non_query_domain_res) %in% query_domains,
                                                                  eval(nodeX_non_query_domain_res)]
        non_query_ids_set1 = all_ids_set1[all_ids_set1 %in% which_set1_proteins]
      } else {
        which_set1_proteins = non_query_domain_res$data_with_pval[non_query_which_nodeY,
                                                                  eval(nodeX_non_query_domain_res)]
        non_query_ids_set1 = all_ids_set1[all_ids_set1 %in% which_set1_proteins]
      }

      # if proteins from domain_res "query" set should not be removed - find those matching domain filtering criteria
      if(!non_query_set_only){
        # if non-query domains should be the same as query domains - do filter by this criteria
        if(query_domains_only){
          query_ids_set1 = all_ids_set1[domain_res$data_with_pval[which_nodeY & eval(nodeZ) %in% query_domains,
                                                                  all_ids_set1 %in% eval(nodeX)]]
        } else {
          query_ids_set1 = all_ids_set1[domain_res$data_with_pval[which_nodeY,
                                                                  all_ids_set1 %in% eval(nodeX)]]
        }
        ids_set1 = unique(c(non_query_ids_set1, query_ids_set1))
      } else ids_set1 = non_query_ids_set1
      # update which datasets should be kept
      to_keep = to_keep &
        length(ids_set1) >= 1 &
        length(non_query_ids_set1) >= non_query_domains_N

      # remove sequences from datasets
      if(to_keep & remove){
        dataset_name = names(interaction_subset[i])
        seq2keep_1 = c(ids_set1, name[2])
        # remove sequences from FASTA sequence list
        fasta_temp = InteractionSubsetFASTA_list$fasta_subset_list[[dataset_name]][seq2keep_1]
        fasta_temp = AAStringSetList(fasta_temp)
        names(fasta_temp) = dataset_name
        InteractionSubsetFASTA_list$fasta_subset_list[dataset_name] = fasta_temp
        # remove interactions from interaction list and update object summaries
        int_temp = InteractionSubsetFASTA_list$interaction_subset[[dataset_name]]
        int_temp$combined_MITAB = int_temp$combined_MITAB[IDs_interactor_B %in% seq2keep_1]
        int_temp$ids_all = seq2keep_1
        int_temp$ids_set1 = int_temp$ids_set1[int_temp$ids_set1 %in% seq2keep_1]
        int_temp$length_set1 = length(int_temp$ids_set1)
        int_temp$length_set2 = length(int_temp$ids_set2)
        InteractionSubsetFASTA_list$interaction_subset[[dataset_name]] = int_temp
      }
    }
    which_to_keep[i] = to_keep
  }

  if(remove){
    # 1 remove datasets
    InteractionSubsetFASTA_list$fasta_subset_list = InteractionSubsetFASTA_list$fasta_subset_list[which_to_keep]
    InteractionSubsetFASTA_list$interaction_subset = InteractionSubsetFASTA_list$interaction_subset[which_to_keep]
    InteractionSubsetFASTA_list$length = length(InteractionSubsetFASTA_list$fasta_subset_list)
    return(InteractionSubsetFASTA_list)
  } else return(which_to_keep)
}
