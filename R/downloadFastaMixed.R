##' Download FASTA for a list of canonical (P04637) and isoform (P04637-2) UniprotAC
##' @name downloadFastaMixed
##' @author Vitalii Kleshchevnikov
##' @param uniprot_ac character vector containg canonical, isoform UniprotAC, or UniprotAC-post-processed chain ID (P04591-PRO_0000261216)
##' @param file_name file name and directory to save the result ("./dir/filename")
##' @return nothing
##' @import Biostrings
##' @export downloadFastaMixed
##' @details \code{downloadFasta} is intended to download from Uniprot and write fasta sequence file (renaming the name of the sequence given by UniProt to UniprotAC only) given generic (P04637) or isoform (P04637-2) UniprotAC. If the file already contains the sequences for given UniprotAC you get a message, if some sequences are missing all sequences will be reloaded.
##' @seealso \code{\link{downloadFasta_postproc}}
##' @examples
##' downloadFastaMixed(uniprot_ac = c("P04637", "P04637-2", "P04591-PRO_0000261216"), file_name = "my_canonical_and_isoform.fasta")
downloadFastaMixed = function(uniprot_ac, file_name){
  # Filtering sequence names by group
  non_canonical = c(grep("-[[:digit:]]+", uniprot_ac, value = F),grep("-PRO_[[:digit:]]+$", uniprot_ac, value = F))
  if(length(non_canonical) >= 1) canonical_all_proteins = uniprot_ac[-non_canonical] else
    canonical_all_proteins = uniprot_ac
  isoform_all_proteins = grep("-[[:digit:]]+", uniprot_ac, value = T)
  canonical_and_isoform = c(canonical_all_proteins, isoform_all_proteins)
  # Filtering sequence names by group
  postproc_all_proteins = grep("-PRO_[[:digit:]]+$", uniprot_ac, value = T)
  postproc = gsub("^[[:alnum:]]+-","",postproc_all_proteins)


  # dowloading FASTA for canonical_and_isoform sequences to a temporary directory
  downloadFasta(uniprot_ac = canonical_and_isoform, file_name = paste0(tempdir(), "/canonical_and_isoform.fasta"))
  # dowloading FASTA for post-processed sequences if file with this data doesn't exist or doesn't contain all sequences
  downloadFastaPostproc(postproc_id = postproc, file_name = paste0(tempdir(), "/postproc.fasta"))

  all_human_viral_proteins.fasta = append(readAAStringSet(paste0(tempdir(), "/canonical_and_isoform.fasta")),
                                          readAAStringSet(paste0(tempdir(), "/postproc.fasta")))
  writeXStringSet(all_human_viral_proteins.fasta,
                  file = file_name,
                  format="fasta")
  return(all_human_viral_proteins.fasta)
}
