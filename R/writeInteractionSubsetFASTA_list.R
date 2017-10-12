##' Write InteractionSubsetFASTA_list into directory system (creating directories as necessary) and return file locations for QSLIMFinder
##' @name writeInteractionSubsetFASTA_list
##' @author Vitalii Kleshchevnikov
##' @param interactionFASTA_list object of class InteractionSubsetFASTA_list created by \code{\link{listInteractionSubsetFASTA}}, and (potentially) filtered using \code{\link{filterInteractionSubsetFASTA_list}} or \code{\link{domainProteinPairMatch}}
##' @param dir directory where to write interactionFASTA_list
##' @return data.table containing path to files and directories for QSLIMFinder: fastafile, queryfile, outputdir, outputfile
##' @import data.table
##' @import Biostrings
##' @export writeInteractionSubsetFASTA_list
##' @seealso \code{\link{listInteractionSubsetFASTA}}
##' @examples
##' forSLIMFinder_file_list = writeInteractionSubsetFASTA_list(interactionFASTA_list = forSLIMFinder_Ready,
##'                            dir = "./SLIMFinder/")
writeInteractionSubsetFASTA_list = function(interactionFASTA_list, dir = "./SLIMFinder/"){
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

  for (i in 1:length(interactionFASTA_list$fasta_subset_list)) {
    name = names(interactionFASTA_list$fasta_subset_list[i])
    name = gsub("\\:", "\\.", name)
    name2 = unlist(split(gsub("interactors_of\\.|\\.","", name), "\\:"))
    fastafile = paste0(input_fasta_dir, name,"fas")
    queryfile = paste0(input_query_dir, name,"fas")
    outputdir = paste0(output_dir, name, "/")
    outputfile = paste0(outputdir, "main_result")

    temp_list = data.table(fastafile = fastafile, queryfile = queryfile,
                           outputdir = outputdir, outputfile = outputfile,
                           interactors_of = name2[1], QSLIMFinder_query = name2[2],
                           sequences = paste0(names(interactionFASTA_list$fasta_subset_list[[i]]),
                                              collapse = "."))
    file_list = rbind(file_list, temp_list)

    writeXStringSet(interactionFASTA_list$fasta_subset_list[[i]],
                    file = fastafile,
                    format="fasta")
    fwrite(list(interactionFASTA_list$interaction_subset[[i]]$ids_set2), queryfile)
  }
  return(file_list)
}
