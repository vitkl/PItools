##' Read and write QSLIMFinder results
##' @rdname readQSLIMFinderMain
##' @name readQSLIMFinderMain
##' @author Vitalii Kleshchevnikov
##' @param outputfile character vector of paths to QSLIMFinder main output files, one way to get these is \code{\link{writeInteractionSubsetFASTA_list}}()$outputfile
##' @param check_exist if TRUE, check that each file exists / is not empty, if false - read all existing and non-empty files
##' @return readQSLIMFinderMain: data.table combined from all QSLIMFinder main output files. Contains summary statistic on each motif found in each dataset
##' @import data.table
##' @export readQSLIMFinderMain
##' @seealso \code{\link{writeInteractionSubsetFASTA_list}}, \code{\link{runQSLIMFinder}}
readQSLIMFinderMain = function(outputfile = "forSLIMFinder_file_list$outputfile", check_exist = F){
  main_result = lapply(outputfile, function(file) {
    #fread("./SLIMFinder/output/interactors_of.O60506_P03496./main_result", stringsAsFactors = F)
    if(check_exist){
      if(file.exists(file)) {
        if(file.size(file) > 0) exist = TRUE else exist = FALSE
      } else exist = FALSE
    } else if(file.exists(file)) if(file.size(file) > 0) fread(file, stringsAsFactors = F)
  })
  if(mean(sapply(main_result, is.null)) == 1) stop("no main_result files found (probably means mis-specified directory)")
  Reduce(rbind, main_result)
}

##' Read and write QSLIMFinder results
##' @rdname readQSLIMFinderMain
##' @name readQSLIMFinderOccurence
##' @author Vitalii Kleshchevnikov
##' @param outputdir character vector of paths to QSLIMFinder output directories, one way to get these is \code{\link{writeInteractionSubsetFASTA_list}}()$outputdir
##' @return readQSLIMFinderOccurence: data.table combined from all QSLIMFinder occurence output files. Contains positions of motifs in each protein.
##' @import data.table
##' @export readQSLIMFinderOccurence
readQSLIMFinderOccurence = function(outputdir = "forSLIMFinder_file_list$outputdir", check_exist = F){
  Occurence_list = lapply(outputdir, function(dir) {
    if(dir.exists(dir)){
      all_files = list.files(dir)
      occurence_file = grep(".occ.csv", all_files, value = T)
      file = paste0(dir,occurence_file)
      if(check_exist){
        if(file.exists(file)) {
          if(file.size(file) > 0) exist = TRUE else exist = FALSE
        } else exist = FALSE
      } else if(length(occurence_file) == 1) if(file.size(file) > 0) fread(file, header = T, stringsAsFactors = F)
    }
  })
  if(mean(sapply(Occurence_list, is.null)) == 1) {
    message(paste0("no occurence files (.occ.csv) found in ",outputdir," (probably means that no motifs were found or mis-specified directory)"))
    return(NULL)
    } else {
    Occurence_list = Occurence_list[!sapply(Occurence_list, is.null)]
    Occurence_list = Occurence_list[sapply(Occurence_list, ncol) == 14]
    # correcting colnames in Occurence_list in case some files have columns with no names or named differently
    Occurence_list = lapply(Occurence_list, function(occurence, Occurence_list_1){
      colnames(occurence) = colnames(Occurence_list_1)
      occurence
    }, Occurence_list[[1]])
    return(Reduce(rbind, Occurence_list))
  }
}

##' Read and write QSLIMFinder results
##' @rdname readQSLIMFinderMain
##' @name writePatternList
##' @author Vitalii Kleshchevnikov
##' @param QSLIMFinder_main_result data.table, the output of readQSLIMFinderMain()
##' @param filename character, paths to file where to write motif list
##' @return writePatternList: data.table containing Name and Pattern columns is saved as a tsv to \code{filename}, a format required for by comparimotif
##' @import data.table
##' @export writePatternList
writePatternList = function(QSLIMFinder_main_result, filename = "./motifs.txt") {
  motifs = unique(QSLIMFinder_main_result[,.(Dataset,Pattern)])
  motifs[, Name := paste0(Dataset,Pattern)][,Dataset:=NULL]
  fwrite(motifs, filename, sep = "\t")
  return(motifs)
}
