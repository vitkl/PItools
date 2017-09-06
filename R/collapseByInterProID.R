##' Collapse InterProScan output to one protein sequence feature per InterProID (or other feature attribute)
##' @name collapseByInterProID
##' @author Vitalii Kleshchevnikov
##' @param InterProScan_features GRanges, the output of \code{readInterProGFF3}, \code{SubsetByInterProEntryType} or \code{addInterProEntryTypes} function
##' @param id_col character (1L) containing the name of the column which defines a feature attribute per which to collapse features (combination of Protein ID and this feature is used)
##' @return \link[GenomicRanges]{GRanges-class} object containing InterProScan output with all protein sequence features collapsed per InterProID (or other feature attribute). Ranges are reduced \link[GenomicRanges]{reduce}, metadata is collapsed in one line and separated with "|".
##' @import Biostrings
##' @import rtracklayer
##' @seealso \link[MItools]{readInterProGFF3}
##' @export collapseByInterProID
##' @examples
##' InterProScan_domains_nonred = collapseByInterProID(InterProScan_features = InterProScan_domains, id_col = "Dbxref")
collapseByInterProID = function(InterProScan_features, id_col = "Dbxref"){
  InterProScan_features_list = split(InterProScan_features, paste0(seqnames(InterProScan_features),"|", mcols(InterProScan_features)[,id_col[1]]))
  InterProScan_features_nonred = GRanges()

  for (i in 1:length(InterProScan_features_list)) {
    Nseq = length(InterProScan_features_list[[i]])
    reduced = reduce(InterProScan_features_list[[i]], with.revmap=T)
    metadata_orig = mcols(InterProScan_features_list[[i]])

    cols2remove = which(colnames(metadata_orig) %in% c("status", "signature_desc", "md5", "ID", "type", "phase"))
    if(length(cols2remove) > 0) metadata_new = metadata_orig[-cols2remove] else metadata_new = metadata_orig
    metadata_new$source = as.character(metadata_new$source)

    mapping = character(Nseq)
    for (ind in 1:length(reduced$revmap)) {
      number = reduced$revmap[[ind]]
      mapping[number] = ind
    }
    metadata_new_list = split(metadata_new, mapping)
    new_metadata_new_list = SplitDataFrameList()
    for (indic in 1:length(metadata_new_list)) {
      if(length(unique(metadata_new$ENTRY_TYPE)) > 1) stop("trying to collapse InterPro sequence features of multiple types! How does it makes sense?")
      temp_metadata = metadata_new_list[[indic]]
      temp_metadata$source = paste0(temp_metadata$source, collapse="|")
      if(id_col != "Dbxref") temp_metadata$Dbxref = paste0(temp_metadata$Dbxref, collapse="|")
      temp_metadata$score = paste0(temp_metadata$score, collapse="|")
      temp_metadata$Target = paste0(temp_metadata$Target, collapse="|")
      temp_metadata$Name = paste0(temp_metadata$Name, collapse="|")
      # convert CharacterList to character first
      temp_metadata$Ontology_term = sapply(temp_metadata$Ontology_term, paste0, collapse=",")
      temp_metadata$Ontology_term = paste0(temp_metadata$Ontology_term, collapse="|")
      new_metadata_new_list[[indic]] = unique(temp_metadata)
    }
    metadata_new = do.call("rbind", new_metadata_new_list)
    mcols(reduced) = metadata_new
    InterProScan_features_nonred = c(InterProScan_features_nonred, reduced)
  }
  # name the vector
  names(InterProScan_features_nonred) = seqnames(InterProScan_features_nonred)
  unique(InterProScan_features_nonred)
}
