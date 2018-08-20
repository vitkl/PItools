##' Read QSLIMFinder (SLIMFinder) occurence output to Genomic Ranges object
##' @rdname SLIMFinderOcc2GRanges
##' @name SLIMFinderOcc2GRanges
##' @author Vitalii Kleshchevnikov
##' @param occurence_file a path to a tsv (txt) file containing QSLIMFinder (SLIMFinder) occurence output
##' @param main_file a path to a tsv (txt) file containing QSLIMFinder (SLIMFinder) main output
##' @param one_from_cloud pick only one motif per motif cloud by the lowest p-value (before multiple hypothesis testing correction)
##' @return Genomic Ranges object containing QSLIMFinder (SLIMFinder) occurence output
##' @import data.table
##' @importFrom GenomicRanges makeGRangesFromDataFrame
##' @export SLIMFinderOcc2GRanges
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{GRangesINinteractionSubsetFASTA}}
SLIMFinderOcc2GRanges = function(occurence_file = "../viral_project/SLIMFinder_Vidal/result/occurence.txt", main_file = "../viral_project/SLIMFinder_Vidal/result/main_result.txt", one_from_cloud = T) {
  occurence = unique(fread(occurence_file, stringsAsFactors = F))
  occurence[, Seq := gsub("_UNK__.+$","",Seq)]
  occurence[, prot_names := Seq]
  if(!is.null(main_file)) {
    pattens = unique(fread(main_file, stringsAsFactors = F)[, c("RunID", "RunTime") := NULL])
    if(one_from_cloud){
      pattens[, order_in_cloud := order(Sig, decreasing = F), by = .(Dataset, Cloud)]
      pattens = unique(pattens[order_in_cloud == 1][, c("order_in_cloud") := NULL])
      pattens = pattens[, Sig_FDR := p.adjust(Sig, "fdr")]
      occurence = occurence[pattens, on = c("Dataset", "Pattern", "Rank", "Sig"), nomatch = 0]
    }
  }

  occurence[, query := gsub("interactors_of\\.[[:alnum:]]{6,10}\\.", "", Dataset)]
  occurence[, interacts_with := gsub("interactors_of\\.|((\\.[[:alnum:]]{6,10})|(\\.[[:alnum:]]{6,10})-PRO_[[:alnum:]]{6,10})", "", Dataset)]
  seqinfo_temp = unique(occurence[,.(Seq, Prot_Len)])
  seqinfo_res = seqinfo_temp$Prot_Len
  names(seqinfo_res) = seqinfo_temp$Seq

  GenomicRanges::makeGRangesFromDataFrame(occurence,
                           keep.extra.columns=T,
                           ignore.strand=T,
                           seqinfo=seqinfo_res,
                           seqnames.field="Seq",
                           start.field="Start_Pos",
                           end.field="End_Pos",
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)
}

##' Merge data.table to GRanges metadata
##' @rdname merge2GRangesmcols
##' @name merge2GRangesmcols
##' @author Vitalii Kleshchevnikov
##' @param range a path to a tsv (txt) file containing QSLIMFinder (SLIMFinder) occurence output
##' @param datatable a path to a tsv (txt) file containing QSLIMFinder (SLIMFinder) main output
##' @param one_from_cloud pick only one motif per motif cloud by the lowest p-value (before multiple hypothesis testing correction)
##' @return Genomic Ranges object containing QSLIMFinder (SLIMFinder) occurence output
##' @import data.table
##' @importFrom GenomicRanges makeGRangesFromDataFrame
##' @importFrom GenomicRanges merge
##' @export merge2GRangesmcols
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{SLIMFinderOcc2GRanges}}
merge2GRangesmcols = function(range, datatable,
                              by.x = c("query", "interacts_with"),
                              by.y = c("IDs_interactor_viral", "IDs_interactor_human")){
  occurence_query_with_domains = GenomicRanges::merge(x = range,
                                       y = datatable,
                                       by.x = by.x,
                                       by.y = by.y,
                                       all.x = T, all.y = F)
  occurence_query_with_domains = unique(occurence_query_with_domains)

  occurence_query_with_domains$width = NULL
  occurence_query_with_domains$strand = NULL
  occurence_query_with_domains = GenomicRanges::makeGRangesFromDataFrame(occurence_query_with_domains,
                                                                         keep.extra.columns=T,
                                                                         ignore.strand=T,
                                                                         seqinfo=seqinfo(range),
                                                                         seqnames.field="seqnames",
                                                                         start.field="start",
                                                                         end.field="end",
                                                                         starts.in.df.are.0based=FALSE)
}

##' Download and read ELM database motif occurrences
##' @rdname ELMdb2GRanges
##' @name ELMdb2GRanges
##' @author Vitalii Kleshchevnikov
##' @param dbfile a path to a gff (txt) file containing ELM database motif occurrences
##' @param dburl url where to get ELM database containing motif occurrences
##' @details ELM database containing motif occurrences in human proteins (9606): "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic="
##' @details ELM database containing motif occurrences in viral proteins (10239): "http://elm.eu.org/instances.gff?q=all&taxon=irus&instance_logic="
##' @return Genomic Ranges object containing ELM database motif occurrences
##' @import GenomicRanges
##' @import data.table
##' @importFrom rtracklayer import.gff3
##' @export ELMdb2GRanges
##' @seealso \code{\link{SLIMFinderOcc2GRanges}}, \code{\link{GRangesINinteractionSubsetFASTA}}
##' @examples
##' instances9606 = ELMdb2GRanges(dbfile = "../viral_project/data_files/instances9606.gff",
##'                   dburl = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic=")
##' instances10239 = ELMdb2GRanges(dbfile = "../viral_project/data_files/instances10239.gff",
##'                   dburl = "http://elm.eu.org/instances.gff?q=all&taxon=irus&instance_logic=")
ELMdb2GRanges = function(dbfile = "../viral_project/data_files/instances9606.gff", dburl = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic=", tsvurl =  "http://elm.eu.org/instances.tsv?q=None&taxon=Homo%20sapiens&instance_logic=", tsvfile = "../viral_project/data_files/instances9606.tsv") {
  if(!file.exists(dbfile)) download.file(dburl, dbfile)
  instances = rtracklayer::import.gff3(dbfile)
  instances$for_benchmarking = 1
  instances$score = NULL
  instances$phase = NULL
  if(!is.null(tsvurl) & !is.null(tsvfile)) {
    if(!file.exists(tsvfile)) download.file(tsvurl, tsvfile)
    instances_tsv = as.data.table(read.delim(tsvfile, skip = 5, stringsAsFactors = F))
    granges_id = paste0(seqnames(instances), start(instances), instances$ID)
    tsv_id = paste0(instances_tsv$Primary_Acc, instances_tsv$Start, instances_tsv$ELMIdentifier)
    instances_tsv = instances_tsv[match(granges_id, tsv_id),]
    mcols(instances) = cbind(mcols(instances), instances_tsv[,.(Accession, ELMType, ProteinName, Primary_Acc, Accessions, References, Methods, InstanceLogic, PDB, Organism)])
  }
  instances
}

##' Filter ELM database motif occurrence GRanges by motif type
##' @rdname filterBYmotifType
##' @name filterBYmotifType
##' @author Vitalii Kleshchevnikov
##' @param instances ELM database motif occurrence GRanges object returned by \code{\link{ELMdb2GRanges}}
##' @param motif_types character vector of motif types
##' @return Genomic Ranges object filtered by motif type (ELMType column)
##' @import GenomicRanges
##' @export filterBYmotifType
##' @seealso \code{\link{SLIMFinderOcc2GRanges}}, \code{\link{GRangesINinteractionSubsetFASTA}}
##' @examples
##' instances9606 = ELMdb2GRanges(dbfile = "../viral_project/data_files/instances9606.gff",
##'                   dburl = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic=")
##' instances9606 = filterBYmotifType(instances9606, motif_types = c("DOC", "MOD", "LIG", "DEG", "CLV", "TRG"))
filterBYmotifType = function(instances, motif_types = c("DOC", "MOD", "LIG", "DEG", "CLV", "TRG")) {
  #TOkeep = lapply(motif_types, function(motif_type) {
  #  grepl(motif_type, instances$ID)
  #})
  #TOkeep = Reduce(`|`,TOkeep)
  TOkeep = instances$ELMType %in% motif_types
  instances[TOkeep]
}
