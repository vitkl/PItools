##' Download FASTA for a list of canonical (P04637) and isoform (P04637-2) UniprotAC
##' @name downloadFasta
##' @author Vitalii Kleshchevnikov
##' @param uniprot_ac character vector containg canonical and isoform UniprotAC
##' @param file_name file name and directory to save the result ("./dir/filename")
##' @return nothing
##' @import Biostrings
##' @export downloadFasta
##' @details \code{downloadFasta} is intended to download from Uniprot and write fasta sequence file (renaming the name of the sequence given by UniProt to UniprotAC only) given generic (P04637) or isoform (P04637-2) UniprotAC. If the file already contains the sequences for given UniprotAC you get a message, if some sequences are missing all sequences will be reloaded.
##' @seealso \code{\link{downloadFasta_postproc}}
##' @examples
##' downloadFasta(uniprot_ac = c("P04637", "P04637-2"), file_name = "my_canonical_and_isoform.fasta")
downloadFasta = function(uniprot_ac, file_name){
  if(if(file.exists(file_name)) mean(uniprot_ac %in% fasta.index(file_name)$desc) == 1 else F){
    message("all sequences for given uniprot_ac are alredy downloaded")
  }
  if(if(file.exists(file_name)) mean(uniprot_ac %in% fasta.index(file_name)$desc) != 1 else T){
    tempdirectory = tempdir()
    message("downloading fasta sequences for given uniprot_ac ...")

    # download all fasta in big files
    sp = loadAllFASTA(uniprot_url = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", uniprot_ac, file = "sp.fasta.gz")
    # tr = loadAllFASTA(uniprot_url = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz", uniprot_ac)
    iso = loadAllFASTA(uniprot_url = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz", uniprot_ac)
    all_fasta = append(sp, iso)
    #all_fasta = append(all_fasta, tr)
    uniprot_ac2 = uniprot_ac[!uniprot_ac %in% names(all_fasta)]

    for(protein in uniprot_ac2){
        file = paste0(tempdirectory, protein)
        tryCatch({
        download.file(paste0("http://www.uniprot.org/uniprot/", protein,".fasta"), file, quiet = T)
        new_fasta = readAAStringSet(file)
        }, error = function(e) e)
        unlink(file)
        # if the search has yielded fasta file with 1 sequence then rename to UniprotAC, otherwise save original name and deal with the problem manually
        names(new_fasta) = gsub("^((sp)|(tr))\\||\\|.+$","",names(new_fasta))
        all_fasta = append(all_fasta, new_fasta)

        if(which(uniprot_ac == protein) %in% seq(1,length(uniprot_ac),10)){
          # write data at each 10th sequence in case R crashes
          writeXStringSet(all_fasta, file = file_name, format="fasta")
          # limit the number of requests per second to 10
          Sys.sleep(1)
        }
    }
    # write final result
    writeXStringSet(all_fasta, file = file_name, format="fasta")
    unlink(tempdirectory)
  }
}

##' Load all FASTA from single UniProt file and filter by uniprot_ac
##' @name downloadFasta
##' @author Vitalii Kleshchevnikov
##' @param uniprot_ac character vector containg canonical and isoform UniprotAC
##' @param uniprot_url file name and directory to save the result ("./dir/filename")
##' @import Biostrings
##' @return AAStringSet
loadAllFASTA = function(uniprot_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", uniprot_ac, file = "sp.fasta.gz") {
  tempdirectory = tempdir()
  download.file(uniprot_url, paste0(tempdirectory,file))
  index_sp = fasta.index(paste0(tempdirectory,file))
  index_sp$desc = gsub("^((sp)|(tr))\\||\\|.+$","",index_sp$desc)
  sp = readAAStringSet(index_sp[index_sp$desc %in% uniprot_ac,])
  names(sp) = gsub("^((sp)|(tr))\\||\\|.+$","",names(sp))
  sp
}
