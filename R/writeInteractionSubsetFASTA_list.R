##' Write InteractionSubsetFASTA_list into directory system (creating directories as necessary) and return file locations for QSLIMFinder
##' @name writeInteractionSubsetFASTA_list
##' @author Vitalii Kleshchevnikov
##' @param interactionFASTA_list object of class InteractionSubsetFASTA_list created by \code{\link{listInteractionSubsetFASTA}}, and (potentially) filtered using \code{\link{filterInteractionSubsetFASTA_list}} or \code{\link{domainProteinPairMatch}}
##' @param dir directory where to write interactionFASTA_list
##' @param analysis_type character 1L, qslimfinder or slimfinder. slimfinder doesn't need query: identical datasets with the same query will be removed and query files not written)
##' @return data.table containing path to files and directories for QSLIMFinder: fastafile, queryfile, outputdir, outputfile
##' @import data.table
##' @import Biostrings
##' @export writeInteractionSubsetFASTA_list
##' @seealso \code{\link{listInteractionSubsetFASTA}}
##' @examples
##' forSLIMFinder_file_list = writeInteractionSubsetFASTA_list(interactionFASTA_list = forSLIMFinder_Ready,
##'                            dir = "./SLIMFinder/")
writeInteractionSubsetFASTA_list = function(interactionFASTA_list, dir = "./SLIMFinder/", analysis_type = c("qslimfinder", "slimfinder")[1]){
  if(!analysis_type %in% c("qslimfinder", "slimfinder")) stop("analysis_type should be one of the \"qslimfinder\", \"slimfinder\"")

  input_dir = paste0(dir, "input/")
  input_fasta_dir = paste0(input_dir, "fasta/")
  input_query_dir = paste0(input_dir, "query/")
  output_dir = paste0(dir, "output/")

  if(!dir.exists(input_dir)) dir.create(input_dir)
  if(!dir.exists(input_fasta_dir)) dir.create(input_fasta_dir)
  if(!dir.exists(input_query_dir)) dir.create(input_query_dir)
  if(!dir.exists(output_dir)) dir.create(output_dir)

  file_list = data.table(fastafile = character(), queryfile = character(),
                         outputdir = character(), outputfile = character(),
                         interactors_of = character(), QSLIMFinder_query = character(),
                         sequences = character())

  #if(analysis_type == "slimfinder") {
  #  seq_indices = 1:length(interactionFASTA_list$fasta_subset_list)
  #  identical = list(i = numeric(), j = numeric())
  #  for(i in seq_indices){
  #    for(j in seq_indices){
  #      i_in_j = mean(names(interactionFASTA_list$fasta_subset_list[[i]]) %in% names(interactionFASTA_list$fasta_subset_list[[j]])) == 1
  #      j_in_i = mean(names(interactionFASTA_list$fasta_subset_list[[j]]) %in% names(interactionFASTA_list$fasta_subset_list[[i]])) == 1
  #      if(i != j & i_in_j & j_in_i){
  #        identical$i = c(identical$i, i)
  #        identical$j = c(identical$j, j)
  #      }
  #    }
  #  }
  #  ind = identical$i < identical$j
  #  identical$i = identical$i[ind]
  #  identical$j = identical$j[ind]
  #  to_keep = seq_indices[!seq_indices %in% c(identical$i, identical$j)]
  #  to_keep = unique(c(to_keep, identical$i[!identical$i %in% identical$j]))
  #  interactionFASTA_list$fasta_subset_list = interactionFASTA_list$fasta_subset_list[to_keep]
  #  interactionFASTA_list$interaction_subset = interactionFASTA_list$interaction_subset[to_keep]
  #}

  # set up progress bar
  pb <- progress::progress_bar$new(
    format = "writing datasets to files [:bar] :current/:total eta: :eta",
    total = length(interactionFASTA_list$fasta_subset_list), clear = FALSE, width= 80, show_after = 0)
  for (i in 1:length(interactionFASTA_list$fasta_subset_list)) {
    pb$tick()
    name = names(interactionFASTA_list$fasta_subset_list[i])
    name = gsub("\\:", "\\.", name)
    name2 = unlist(strsplit(gsub("interactors_of\\.|\\.$", "",name), "\\."))
    if(analysis_type == "slimfinder") {
      name2[2] = NA
      name = paste0("interactors_of.", name2[1], ".")
      }
    fastafile = paste0(input_fasta_dir, name,"fas")
    queryfile = paste0(input_query_dir, name,"fas")
    outputdir = paste0(output_dir, name, "/")
    outputfile = paste0(outputdir, "main_result")

    seqnames_ = names(interactionFASTA_list$fasta_subset_list[[i]])
    seqnames_ = seqnames_[order(seqnames_)]

    temp_list = data.table(fastafile = fastafile, queryfile = queryfile,
                           outputdir = outputdir, outputfile = outputfile,
                           interactors_of = name2[1], QSLIMFinder_query = name2[2],
                           sequences = paste0(seqnames_,
                                              collapse = "."))
    file_list = rbind(file_list, temp_list)

    fasta_seq = interactionFASTA_list$fasta_subset_list[[i]]
    if(length(fasta_seq[unique(names(fasta_seq))]) < length(fasta_seq)) {
      fasta_seq = fasta_seq[unique(names(fasta_seq))]
      warning(paste0(temp_list$interactors_of,": duplicated sequences detected and removed (likely reason: set1 and set2 have 1 or more of the same interacting proteins)"))
    }

    writeXStringSet(fasta_seq,
                    file = fastafile,
                    format="fasta")
    if(analysis_type == "qslimfinder") {
      fwrite(list(interactionFASTA_list$interaction_subset[[i]]$ids_set2), queryfile)
    }
  }

  file_list = unique(file_list)
  return(file_list)
}
