##' Download FASTA for a list of post-processed chain ID (e.g. PRO_0000261216)
##' @name download_fasta_postproc
##' @author Vitalii Kleshchevnikov
##' @param uniprot_ac character vector containg post-processed chain ID
##' @param file_name file name and directory to save the result ("./dir/filename")
##' @return nothing
##' @import Biostrings
##' @importFrom GenomeInfoDb seqnames
##' @importMethodsFrom rtracklayer import
##' @export download_fasta_postproc
##' @details function to download from Uniprot and write fasta sequence files (renaming the UniProt-given name of the sequence to proteinID-FeatureID: P04591-PRO_0000261216) given post-processed chain ID, such as PRO_0000261216. If the file already contains the sequences for given post-processed chain ID you get a message, if some sequences are missing all sequences will be reloaded.
##' @seealso \code{\link{download_fasta}}
##' @examples
##' download_fasta_postproc(postproc_id = c("PRO_0000261216"), file_name = "my_post_processed_chain.fasta")
download_fasta_postproc = function(postproc_id, file_name){
  if(if(file.exists(file_name)) mean(postproc_id %in% gsub("^[[:alnum:]]+-","",fasta.index(file_name)$desc)) == 1 else F){
    message("all sequences for given post-processed chain ID are alredy downloaded")
  }
  if(if(file.exists(file_name)) mean(postproc_id %in% gsub("^[[:alnum:]]+-","",fasta.index(file_name)$desc)) != 1 else T){
    message("downloading fasta sequences for given post-processed chain IDs ...")
    postproc_fasta = AAStringSet()
    for(feature in postproc_id){
      # Read gff from UniProt
      new_feature = import(paste0("http://www.uniprot.org/uniprot/?query=",feature,"&format=gff"), format="gff")
      # subset features with ID
      new_feature = new_feature[!is.na(new_feature$ID)]
      # select the feature we have searched for
      new_feature = new_feature[new_feature$ID == feature]
      # read FASTA for the protein
      protein = as.character(seqnames(new_feature))
      proteinFASTA = readAAStringSet(paste0("http://www.uniprot.org/uniprot/", protein,".fasta"))
      # select the first and the only sequence in a set and subset it by feature range
      postproc_fasta_new = AAStringSet(proteinFASTA[[1]][ranges(new_feature)])
      names(postproc_fasta_new) = paste0(protein, "-", feature)
      postproc_fasta = append(postproc_fasta, postproc_fasta_new)

      if(which(postproc_id == feature) %in% seq(1,length(postproc_id),10)){
        # write data at each 10th sequence in case R crashes
        writeXStringSet(postproc_fasta, file = file_name, format="fasta")
        # limit the number of requests per second to 10
        Sys.sleep(1)
      }
    }
    # write final result
    writeXStringSet(postproc_fasta, file = file_name, format="fasta")
  }
}
