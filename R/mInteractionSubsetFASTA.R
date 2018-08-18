##' get FASTA sequences of proteins that interact with seed proteins (datasets for QSLIMFinder)
##' @name listInteractionSubsetFASTA
##' @author Vitalii Kleshchevnikov
##' @param interaction_set1 cleanMItab object containing molecular interactions (such as all interactions between human proteins)
##' @param interaction_set2 a different cleanMItab object containing molecular interactions (conceptually different source of interaction data, such as other type of interactions or human-viral interactions)
##' @param seed_id_vect a vector of interactor IDs for which to retrieve interactions (for each separately)
##' @param fasta AAStringSet containing sequences for all proteins in interaction_set1 and interaction_set2
##' @param single_interact_from_set2 logical, split sequence sets to contain only one protein from interaction_set2 (only one query protein for QSLIMFinder). If FALSE, set2 will contain all proteins that interact with an element of \code{seed_id_vect} (which means multiple query proteins for QSLIMFinder).
##' @param set1_only logical, only relevant if \code{single_interact_from_set2 = TRUE}, sequence set1 should contain only proteins that interact with an element of seed_id_vect in interaction_set1. If FALSE, proteins that interact with an element of seed_id_vect in interaction_set2 but are not a single query protein are also included. Argument for \code{\link{singleInteractFromSet2}}
##' @return object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @import data.table
##' @import Biostrings
##' @export listInteractionSubsetFASTA
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{recodeFASTA}}, \code{\link{interactionSubsetMapID}}, \code{\link{recodeANDinteractionSubsetFASTA}}, \code{\link{singleInteractFromSet2}}, \code{\link{listSingleInteractFromSet2}}, \code{\link{filterInteractionSubsetFASTA_list}}
##' @examples
##' forSLIMFinder = listInteractionSubsetFASTA(interaction_set1 = all_human_interaction,
##'                  interaction_set2 = all_viral_interaction,
##'                  seed_id_vect = proteins_w_signif_domains,
##'                  fasta = all.fasta,
##'                  single_interact_from_set2 = T, set1_only = T)
listInteractionSubsetFASTA = function(interaction_set1, interaction_set2, seed_id_vect,
                                      fasta, single_interact_from_set2 = T, set1_only = T) {

  if(!grepl("clean_MItab",class(interaction_set1))) stop("interaction_set1 is not of class clean_MItab27 or related clean_MItab class")
  if(!grepl("clean_MItab",class(interaction_set2))) stop("interaction_set2 is not of class clean_MItab27 or related clean_MItab class")

  # remove seed proteins with no FASTA
  seed_id_vect = seed_id_vect[seed_id_vect %in% names(fasta)]
  # remove seed proteins with no interactions in both sets
  seed_id_vect = seed_id_vect[seed_id_vect %in% extractInteractors(interaction_set1) &
                                seed_id_vect %in% extractInteractors(interaction_set2)]

  subset1 = subset2setsBy1ID(interaction_set1 = interaction_set1,
                             interaction_set2 = interaction_set2,
                             seed_id = seed_id_vect[1])
  subset1_fasta = listSingleInteractFromSet2(subset1, single_interact_from_set2, set1_only, fasta)

  if(length(seed_id_vect) >= 2){
    # set up progress bar
    pb = progress::progress_bar$new(
      format = "Creating datasets for QSLIMFinder [:bar] :current/:total eta: :eta",
      total = length(seed_id_vect)-1, clear = FALSE, width= 80, show_after = 0)
    for (seed_id in seed_id_vect[2:length(seed_id_vect)]) {
      pb$tick()
      subset = subset2setsBy1ID(interaction_set1 = interaction_set1,
                                interaction_set2 = interaction_set2,
                                seed_id = seed_id)
      subset_fasta = listSingleInteractFromSet2(subset1 = subset, single_interact_from_set2, set1_only, fasta)
      subset1_fasta$fasta_subset_list = c(subset1_fasta$fasta_subset_list, subset_fasta$fasta_subset_list)
      subset1_fasta$interaction_subset = c(subset1_fasta$interaction_subset, subset_fasta$interaction_subset)
    }
  }

  subset1_fasta$length = length(subset1_fasta$fasta_subset_list)
  class(subset1_fasta) = "InteractionSubsetFASTA_list"
  return(subset1_fasta)
}

##' Recode names and retrieve FASTA sequences of interactors of a specific protein
##' @name recodeANDinteractionSubsetFASTA
##' @author Vitalii Kleshchevnikov
##' @param subset1 interaction_subset_from_2_sets object produced by \code{\link{subset2setsBy1ID}}
##' @param fasta AAStringSet, names are UniProtKB accessions
##' @return object of class interactionSubsetFASTA containing: FASTA sequences for interacting proteins, original interaction_subset_from_2_sets
##' @import data.table
##' @import Biostrings
##' @export recodeANDinteractionSubsetFASTA
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{recodeFASTA}}, \code{\link{listInteractionSubsetFASTA}}
##' @examples
##' subset1_fasta = recodeANDinteractionSubsetFASTA(subset1, fasta)
recodeANDinteractionSubsetFASTA = function(subset1, fasta){
  #subset1 = interactionSubsetMapID(subset1, fasta$names_mapping)
  #subset1_fasta = interactionSubsetFASTA(int_subset = subset1, fasta = fasta$fasta)
  subset1_fasta = interactionSubsetFASTA(int_subset = subset1, fasta = fasta)
  #subset1_fasta = recodeFASTA(subset1_fasta, fasta$names_mapping, i = 1)
  return(subset1_fasta)
}

##' Filter interaction subset from 2 sets: keep only i-th element from the set2
##' @name singleInteractFromSet2
##' @author Vitalii Kleshchevnikov
##' @param subset1 interaction_subset_from_2_sets object produced by \code{\link{subset2setsBy1ID}}
##' @param set1_only logical, only relevant if \code{single_interact_from_set2 = TRUE}, sequence set1 should contain only proteins that interact with an element of seed_id_vect in interaction_set1. If FALSE, proteins that interact with an element of seed_id_vect in interaction_set2 with the exception of the single one () are also included
##' @param i which protein in the set2 inside \code{subset1} to keep as a set2 element
##' @return object of class interaction_subset_from_2_sets, the same as returned by \code{\link{subset2setsBy1ID}}, but set2 contains only one protein and (optionally) all other set2 proteins and corresponding interactions are removed
##' @import data.table
##' @import Biostrings
##' @export singleInteractFromSet2
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{recodeFASTA}}, \code{\link{listInteractionSubsetFASTA}}
##' @examples
##' subset1_1 = singleInteractFromSet2(subset1, set1_only = set1_only, i = 1)
singleInteractFromSet2 = function(subset1, set1_only, i){
  subset1_1 = subset1
  subset1_1$name = paste0(subset1$name,":",subset1$ids_set2[i])
  subset1_1$ids_set2 = subset1$ids_set2[i]
  subset1_1$length_set2 = 1
  if(!set1_only){
    if("NULL" %in% subset1_1$ids_set1){
      if("NULL" %in% subset1_1$ids_set2){
        stop("no interactions in both set 1 and set 2")
      } else {
        if(length(subset1$ids_set2[-i]) >= 1){
          subset1_1$ids_set1 = subset1$ids_set2[-i]
          subset1_1$length_set1 = length(subset1_1$ids_set1)
        }
      }
    } else {
      if("NULL" %in% subset1_1$ids_set2){
        subset1_1$ids_all = subset1_1$ids_set1
      } else {
        if(length(subset1$ids_set2[-i]) >= 1){
          subset1_1$ids_set1 = c(subset1$ids_set1, subset1$ids_set2[-i])
          subset1_1$length_set1 = length(subset1_1$ids_set1)
        }
      }
    }
  }
  if("NULL" %in% subset1_1$ids_set1){
    if("NULL" %in% subset1_1$ids_set2){
      stop("no interactions in both set 1 and set 2")
    } else {
      subset1_1$ids_all = subset1_1$ids_set2
    }
  } else {
    if("NULL" %in% subset1_1$ids_set2){
      subset1_1$ids_all = subset1_1$ids_set1
    } else {
      subset1_1$ids_all = c(subset1_1$ids_set1, subset1_1$ids_set2)
    }
  }


  subset1_1$combined_MITAB = subset1_1$combined_MITAB[(IDs_interactor_A == subset1$name &
                                                         IDs_interactor_B %in% subset1_1$ids_all) |
                                                        (IDs_interactor_B == subset1$name &
                                                           IDs_interactor_A %in% subset1_1$ids_all),]
  return(subset1_1)
}

##' Split interaction subset from 2 sets so that each elements contains only one of the proteins from set 2, and retrieve FASTA sequences for each
##' @name listSingleInteractFromSet2
##' @author Vitalii Kleshchevnikov
##' @param subset1 interaction_subset_from_2_sets object produced by \code{\link{subset2setsBy1ID}}
##' @param single_interact_from_set2 logical, split sequence sets to contain only one protein from interaction_set2 (only one query protein for QSLIMFinder). If FALSE, set2 will contain all proteins that interact with an element of \code{seed_id_vect} (which means multiple query proteins for QSLIMFinder).
##' @param set1_only logical, only relevant if \code{single_interact_from_set2 = TRUE}, sequence set1 should contain only proteins that interact with an element of seed_id_vect in interaction_set1. If FALSE, proteins that interact with an element of seed_id_vect in interaction_set2 with the exception of the single one () are also included. Argument for \code{\link{singleInteractFromSet2}}
##' @param fasta list containing AAStringSet renamed by \code{\link{recodeFASTA}} and mapping table from \code{\link{recodeFASTA}}
##' @return object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @import data.table
##' @import Biostrings
##' @export listSingleInteractFromSet2
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{listInteractionSubsetFASTA}}
##' @examples
##' subset1_fasta = listSingleInteractFromSet2(subset1, single_interact_from_set2, set1_only, fasta)
listSingleInteractFromSet2 = function(subset1, single_interact_from_set2, set1_only, fasta){
  if(single_interact_from_set2){
    if(subset1$length_set2 >= 1){
      subset1_1 = singleInteractFromSet2(subset1, set1_only = set1_only, i = 1)
      subset1_1_fasta = interactionSubsetFASTA(int_subset = subset1_1, fasta = fasta)
      subset1_1_fasta$fasta_subset_list = subset1_1_fasta$fasta_subset_list[unique(names(subset1_1_fasta$fasta_subset_list))]
      if(subset1$length_set2 >= 2){
        for (ind in 2:subset1$length_set2) {
          subset1_temp = singleInteractFromSet2(subset1, set1_only = set1_only, i = ind)
          subset1_temp_fasta = interactionSubsetFASTA(int_subset = subset1_temp, fasta = fasta)
          subset1_1_fasta$fasta_subset_list = c(subset1_1_fasta$fasta_subset_list,
                                                subset1_temp_fasta$fasta_subset_list[
                                                  unique(names(subset1_temp_fasta$fasta_subset_list))])
          subset1_1_fasta$interaction_subset = c(subset1_1_fasta$interaction_subset,
                                                 subset1_temp_fasta$interaction_subset)
        }
      }
      subset1_fasta = subset1_1_fasta
    } else {
      warning(paste0("set 2 is empty, name:", subset1$name))
      subset1_fasta = interactionSubsetFASTA(int_subset = subset1, fasta = fasta)
    }
  } else {
    subset1_fasta = interactionSubsetFASTA(int_subset = subset1, fasta = fasta)
  }
  return(subset1_fasta)
}
