##' Center QSLIMFinder datasets (InteractionSubsetFASTA_list) at domains
##' @name centerDomains
##' @author Vitalii Kleshchevnikov
##' @param interactionFASTA_list object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @param domain_res domain enrichment results (object of S3 class "XYZinteration_XZEmpiricalPval" (list), the output of permutationPval() call), or data.table with domain_col, query_col and seed_col
##' @param domain_col name of the column in \code{domain_res} containing domain ID, not relevant if \code{domain_res} is XYZinteration_XZEmpiricalPval
##' @param query_col name of the column in \code{domain_res} containing query protein ID, not relevant if \code{domain_res} is XYZinteration_XZEmpiricalPval
##' @param seed_col name of the column in \code{domain_res} containing seed protein ID, not relevant if \code{domain_res} is XYZinteration_XZEmpiricalPval
##' @return object of class InteractionSubsetFASTA_list containing domain-centered datasets
##' @import data.table
##' @import Biostrings
##' @export centerDomains
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{recodeFASTA}}, \code{\link{interactionSubsetMapID}}, \code{\link{recodeANDinteractionSubsetFASTA}}, \code{\link{singleInteractFromSet2}}, \code{\link{listSingleInteractFromSet2}}, \code{\link{filterInteractionSubsetFASTA_list}}
##' @examples
##' forSLIMFinder = listInteractionSubsetFASTA(interaction_set1 = all_human_interaction,
##'                  interaction_set2 = all_viral_interaction,
##'                  seed_id_vect = proteins_w_signif_domains,
##'                  fasta = all.fasta,
##'                  single_interact_from_set2 = T, set1_only = T)
centerDomains = function(interactionFASTA_list, domain_res,
                         domain_col = "IDs_domain_human",
                         query_col = "IDs_interactor_viral",
                         seed_col = "IDs_interactor_human"){
  # get unique dataset ids
  name = names(interactionFASTA_list$fasta_subset_list)
  name2 = strsplit(gsub("interactors_of\\.|\\.$", "",name), "\\:")
  datasets = data.table(seed = sapply(name2, function(name) name[1]),
                        query = sapply(name2, function(name) name[2]))
  # get domain presence and viral interaction data
  if(class(domain_res) == "XYZinteration_XZEmpiricalPval"){
    domains = unique(domain_res$data_with_pval[,c(domain_res$nodes$nodeZ, domain_res$nodes$nodeX, domain_res$nodes$nodeY), with = F])
  } else if(class(domain_res) == "data.table") {
    domains = unique(domain_res[,c(domain_col, query_col, seed_col), with = F])
  } else stop("domain_res is neither XYZinteration_XZEmpiricalPval nor data.table")

  domain_datasets = merge(x = datasets, y = domains,
                          by.x = c("seed", "query"), by.y = c(seed_col, query_col),
                          all.x = F, all.y = F)
  setnames(domain_datasets, domain_col, "domain")

  # generate instructions to construct domain-centered datasets and dataset names
  domain_datasets_list = split(domain_datasets, paste0("interactors_of.",domain_datasets$domain,":",domain_datasets$query,"."))

  # construct domain-centered datasets
  domainInteractionFASTA_list = list(fasta_subset_list = NULL, interaction_subset = NULL, length = integer())
  for(j in 1:length(domain_datasets_list)) {
    domain_dataset = domain_datasets_list[[j]]
    name = paste0(unique(domain_dataset$domain), ":", unique(domain_dataset$query))
    indx = datasets[, which(seed %in% domain_dataset$seed & query == unique(domain_dataset$query))]

    eval(parse(text = paste0("domainInteractionFASTA_list$fasta_subset_list$`",name,"` = Reduce(c, interactionFASTA_list$fasta_subset_list[indx])")))

    interaction_subset_list = interactionFASTA_list$interaction_subset[indx]
    interaction_subset = interaction_subset_list[[1]]

    # combine interaction datasets for multiple seed proteins
    if(length(interaction_subset_list) >= 2){
      for(i in 2:length(interaction_subset_list)){
        interaction_subset2 = interaction_subset_list[[i]]
        interaction_subset$combined_MITAB = unique(rbind(interaction_subset$combined_MITAB, interaction_subset2$combined_MITAB))
        interaction_subset$ids_all = unique(c(interaction_subset$ids_all, interaction_subset2$ids_all))
        interaction_subset$ids_set1 = unique(c(interaction_subset$ids_set1, interaction_subset2$ids_set1))
        interaction_subset$ids_set2 = unique(c(interaction_subset$ids_set2, interaction_subset2$ids_set2))
      }
    }
    # modify metadata
    interaction_subset$name = name
    interaction_subset$length_set1 = length(interaction_subset$ids_set1)
    interaction_subset$length_set2 = 1
    interaction_subset$description = paste0("Interacting partners of seed_id (",unique(domain_dataset$domain),":", paste0(domain_dataset$seed, collapse = ", "),") from interaction_set1 () and interaction_set2 (", unique(domain_dataset$query),").")
    # add domain annotations
    interaction_subset$domain_in = domain_dataset
    class(interaction_subset) = "interaction_subset_from_2_sets"

    eval(parse(text = paste0("domainInteractionFASTA_list$interaction_subset$`",name,"` = interaction_subset")))
  }
  domainInteractionFASTA_list$length = length(domainInteractionFASTA_list$fasta_subset_list)
  class(domainInteractionFASTA_list) = "InteractionSubsetFASTA_list"
  domainInteractionFASTA_list
}
