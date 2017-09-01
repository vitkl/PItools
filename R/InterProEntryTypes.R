##' Download FASTA for a list of canonical (P04637) and isoform (P04637-2) UniprotAC
##' @name readInterProGFF3
##' @author Vitalii Kleshchevnikov
##' @param filename.gz character, file name and directory to read InterProScan output in the gff format ("./dir/filename.gz")
##' @param InterProScan_result GRanges, the output of \code{readInterProGFF3} function
##' @param selected_ENTRY_TYPE character, select InterPro signature matches from \code{InterProScan_result} using these InterPro Entry Types
##' @param entry.list_path character, file name and directory where to store InterPro sca
##' @return \code{readInterProGFF3}, \code{addInterProEntryTypes} or \code{SubsetByInterProEntryType}: \link[GenomicRanges]{GRanges-class} object containing InterProScan output appended by InterPro Entry Types information, with names metadata and sequence length information imported correctly
##' @return \code{getInterProEntryTypes}: data.table containing InterPro Entry List \link{ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list}
##' @return \code{getInterPro2memberDB}: data.table containing  InterProID to member database ID mapping
##' @import Biostrings
##' @import rtracklayer
##' @export readInterProGFF3
##' @export addInterProEntryTypes
##' @export getInterProEntryTypes
##' @export SubsetByInterProEntryType
##' @export getInterPro2memberDB
##' @details For further details on the InterProScan output format please visit the link \link{https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats}
##' @examples
##' addInterProEntryTypes(filename.gz = "./processed_data_files/all_human_viral_protein_domains.gff3.gz")
addInterProEntryTypes = function(InterProScan_result, entry.list_path){
  # 1 add InterPro Entry Types information
  # download and/or read InterPro_entry_types mapping table
  InterPro_entry_types = getInterProEntryTypes(entry.list_path)
  # find matching signatures
  matching_signatures = match(as.character(InterProScan_result$Dbxref), InterPro_entry_types$ENTRY_AC)
  # add signature type information to the gff2 file
  InterProScan_result$ENTRY_TYPE = InterPro_entry_types$ENTRY_TYPE[matching_signatures]
  return(InterProScan_result)
}

readInterProGFF3 = function(filename.gz){
  # 1 ungzip and read InterProScan result
  filename = substr(filename.gz, 1, nchar(filename.gz)-3)
  gunzip(filename.gz, remove = F)
  InterProScan_result = import(con = filename, format = "gff3")
  unlink(filename)
  # clean InterPro ID in the Dbxref column
  InterProScan_result$Dbxref = gsub("InterPro:","", InterProScan_result$Dbxref)
  InterProScan_result$Dbxref = gsub("\"","", InterProScan_result$Dbxref)
  # InterProScan_result$Ontology_term = gsub("\"","", InterProScan_result$Ontology_term)

  # 2 add names metadata to InterProScan_result allow subsetting with a character
  names(InterProScan_result) = seqnames(InterProScan_result)

  # 3 import seqence length information correctly into the GRanges object format
  # find seqence length
  seq_length = end(InterProScan_result)[InterProScan_result$type == "polypeptide"]
  # name the vector
  names(seq_length) = seqnames(InterProScan_result)[InterProScan_result$type == "polypeptide"]
  # match names in the vector and in the GRanges seqlengths slot
  seqlengths(InterProScan_result) = seq_length[match(names(seqlengths(InterProScan_result)),names(seq_length))]
  return(InterProScan_result)
}

SubsetByInterProEntryType = function(InterProScan_result, select_ENTRY_TYPE = c("Domain", "Active_site", "Binding_site", "Conserved_site", "PTM", "Repeat")){
  # create a subset that contains "Domain", "Active_site", "Binding_site", "Conserved_site", "PTM" signatures
  InterProScan_domains = InterProScan_result[InterProScan_result$ENTRY_TYPE %in% select_ENTRY_TYPE]
}

getInterProEntryTypes = function(entry.list_path){
  if(!file.exists(entry.list_path)) download("ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list", entry.list_path)
  InterPro_entry_types = fread(entry.list_path, stringsAsFactors = F)
  return(InterPro_entry_types)
}

getInterPro2memberDB = function(InterProScan_result){
  unique(data.table(InterProID = InterProScan_domains$Dbxref, memberDBID = InterProScan_domains$Name))
}
