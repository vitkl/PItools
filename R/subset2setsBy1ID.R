  ##' Subset 2 sets of protein interaction data using a single ID and combine the results
  ##' @name subset2setsBy1ID
  ##' @author Vitalii Kleshchevnikov
  ##' @param interaction_set1 cleanMItab object containing molecular interactions
  ##' @param interaction_set2 different cleanMItab object containing molecular interactions
  ##' @param seed_id interactor ID for which to retrieve interactions
  ##' @return object of class interaction_subset_from_2_sets (list): name = seed_id,
  ##' @return combined_MITAB = combined interaction data for seed_id from interaction_set1 and interaction_set2
  ##' @return ids_all = IDs of all interactors of seed_id
  ##' @return ids_set1 = IDs of interactors of seed_id in interaction_set1
  ##' @return ids_set2 = IDs of interactors of seed_id in interaction_set1
  ##' @return length_set1, length_set2 - number of interactions in each set
  ##' @return description - short decription of the data
  ##' @import data.table
  ##' @export subset2setsBy1ID
  ##' @export print.interaction_subset_from_2_sets
  ##' @seealso \code{\link{interactionSubsetMapID}}
  ##' @examples
  ##' subset = subset2setsBy1ID(interaction_set1 = all_human_interaction,
  ##'                         interaction_set2 = all_viral_interaction,
  ##'                         seed_id = proteins_w_signif_domains[2])
  ##' # print the result
  ##' subset
subset2setsBy1ID = function(interaction_set1, interaction_set2, seed_id){
  if(length(seed_id) != 1) warning("WARNING: multiple ids may not be meaningful")
  from_set1 = subsetMITABbyID(interaction_set1, ID_seed = seed_id,
                              within_seed = F, only_seed2nonseed = T)
  from_set2 = subsetMITABbyID(interaction_set2, ID_seed = seed_id,
                              within_seed = F, only_seed2nonseed = T)
  length_set1 = length(unique(from_set1$data$pair_id))
  length_set2 = length(unique(from_set2$data$pair_id))

  if(length_set1 > 0) ids_set1 = unique(unlist(from_set1$data[,IDs_interactor_B])) else ids_set1 = "NULL"
  if(length_set2 > 0) ids_set2 = unique(unlist(from_set2$data[,IDs_interactor_B])) else ids_set2 = "NULL"

  if(length_set1 > 0) {
    if(length_set2 > 0) ids_all = c(ids_set1, ids_set2) else ids_all = c(ids_set1)} else {
      if(length_set2 > 0) ids_all = c(ids_set2) else stop(paste0("seed_id (",seed_id,") has no interactions is both interaction_set1 and interaction_set2"))
    }

  combined = rbind(from_set1$data, from_set2$data)
  res = list(name = seed_id, combined_MITAB = combined, ids_all = ids_all,
             ids_set1 = ids_set1, ids_set2 = ids_set2,
             length_set1 = length_set1, length_set2 = length_set2,
             description = paste0("Interacting partners of seed_id (",seed_id,") from interaction_set1 (", match.call()[["interaction_set1"]],") and interaction_set2 (", match.call()[["interaction_set2"]],")."))
  class(res) = "interaction_subset_from_2_sets"
  return(res)
}

print.interaction_subset_from_2_sets = function(data){
  cat("\n\n")
  cat(data$description)
  cat("\n\n This format is useful for preparing data for QSLIMFinder where one set is the search set and the other set is the query set. \n\n")
  cat(paste0("set 1 identifiers (length ",data$length_set1,"): \n"))
  print(data$ids_set1)
  cat(paste0("set 2 identifiers (length ",data$length_set2,"): \n"))
  print(data$ids_set2)
  cat("\n data in MITAB format is at $combined_MITAB (data.table)")
  cat(paste0("\n object class (S3): ",class(data)))
}

##' interaction_subset_from_2_sets renames ids in the interaction_subset_from_2_sets object
##' @name interactionSubsetMapID
##' @author Vitalii Kleshchevnikov
##' @param subset interaction_subset_from_2_sets, the output of \code{\link{subset2setsBy1ID}}
##' @param mapping how to translate original IDs to new IDs
##' @return object of class interaction_subset_from_2_sets (list) with ids_all slot containing new names and slot ids_all_old containing original names
##' @import data.table
##' @export interactionSubsetMapID
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{recodeFASTA}}
##' @examples
##' subset = interactionSubsetMapID(subset, all.fasta$names_mapping)
interactionSubsetMapID = function(subset, mapping) {
  subset$ids_all_old = subset$ids_all
  subset$ids_all = unique(mapping[old_names %in% subset$ids_all_old, new_names])
  return(subset)
}
