##' FASTA sequences of interactors of a specific protein
##' @name interactionSubsetFASTA
##' @author Vitalii Kleshchevnikov
##' @param int_subset interaction_subset_from_2_sets object renamed using interactionSubsetMapID
##' @param fasta AAStringSet renamed by recodeFASTA
##' @return object of class interactionSubsetFASTA containing: FASTA sequences for interacting proteins, original interaction_subset_from_2_sets
##' @import data.table
##' @import Biostrings
##' @export interactionSubsetFASTA
##' @seealso \code{\link{subset2setsBy1ID}}, \code{\link{recodeFASTA}}, \code{\link{interactionSubsetMapID}}
##' @examples
##' subset = subset2setsBy1ID(interaction_set1 = all_human_interaction,
##'                         interaction_set2 = all_viral_interaction,
##'                         seed_id = proteins_w_signif_domains[2])
##' # subset FASTA
##' subset_fasta = interactionSubsetFASTA(int_subset = subset, fasta = all.fasta$fasta)
##' # rename FASTA sequences to original names
##' subset_fasta = recodeFASTA(subset_fasta, all.fasta$names_mapping, i = 1)
interactionSubsetFASTA = function(int_subset, fasta){
  if(class(int_subset) != "interaction_subset_from_2_sets") stop("Invalid input, int_subset should be the output of subset2setsBy1ID(), class interaction_subset_from_2_sets")
  if(class(fasta) != "AAStringSet") stop(" fasta should be of class AAStringSet (Biostrings) ")
  fasta_subset = fasta[int_subset$ids_all]
  fasta_subset_list = AAStringSetList(fasta_subset)
  names(fasta_subset_list) = paste0("interactors_of.",int_subset$name,".")
  int_subset_list = list(int_subset)
  names(int_subset_list) = paste0("interactors_of.",int_subset$name,".")
  res = list(fasta_subset_list = fasta_subset_list, interaction_subset = int_subset_list)
  class(res) = "interactionSubsetFASTA"
  return(res)
}
