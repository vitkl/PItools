##' Download FASTA for a list of canonical (P04637) and isoform (P04637-2) UniprotAC
##' @name download_fasta
##' @author Vitalii Kleshchevnikov
##' @param uniprot_ac character vector containg canonical and isoform UniprotAC
##' @param file_name file name and directory to save the result ("./dir/filename")
##' @return nothing
##' @import Biostrings
##' @export download_fasta
##' @details \code{download_fasta} is intended to download from Uniprot and write fasta sequence file (renaming the name of the sequence given by UniProt to UniprotAC only) given generic (P04637) or isoform (P04637-2) UniprotAC. If the file already contains the sequences for given UniprotAC you get a message, if some sequences are missing all sequences will be reloaded.
##' @seealso \code{\link{download_fasta_postproc}}
##' @examples
##' download_fasta(uniprot_ac = c("P04637", "P04637-2"), file_name = "my_canonical_and_isoform.fasta")
download_fasta = function(uniprot_ac, file_name){
  if(if(file.exists(file_name)) mean(uniprot_ac %in% fasta.index(file_name)$desc) == 1 else F){
    message("all sequences for given uniprot_ac are alredy downloaded")
  }
  if(if(file.exists(file_name)) mean(uniprot_ac %in% fasta.index(file_name)$desc) != 1 else T){
    message("downloading fasta sequences for given uniprot_ac ...")
    all_fasta = AAStringSet()
    sequences_per_id = numeric(length = length(uniprot_ac))
    names(sequences_per_id) = uniprot_ac
    for(protein in uniprot_ac){
      new_fasta = readAAStringSet(paste0("http://www.uniprot.org/uniprot/", protein,".fasta"))
      # if the search has yielded fasta file with 1 sequence then rename to UniprotAC, otherwise save original name and deal with the problem manually
      if(length(new_fasta) == 1) names(new_fasta) = protein
      all_fasta = append(all_fasta, new_fasta)
      sequences_per_id[protein] = length(new_fasta)

      if(which(uniprot_ac == protein) %in% seq(1,length(uniprot_ac),10)){
        # write data at each 10th sequence in case R crashes
        writeXStringSet(all_fasta, file = file_name, format="fasta")
        # limit the number of requests per second to 10
        Sys.sleep(1)
      }
    }
    # write final result
    writeXStringSet(all_fasta, file = file_name, format="fasta")
  }
}
