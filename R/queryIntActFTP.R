##' load all data from IntAct ftp
##' @name loadIntActFTP
##' @author Vitalii Kleshchevnikov
##' @description loadIntActFTP function loads all data that is stored in \link{ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/} in intact.txt file. Details: \link{http://www.ebi.ac.uk/intact/downloads}. This file contains the data from the following databases "IntAct", "MINT", "DIP", "bhf-ucl", "MPIDB", "MatrixDB", "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo"
##' @param file directory where to save/look for the local copy
##' @return saves intact.txt to a file, returns the object of class RAW_MItab27
##' @import data.table
##' @importFrom R.utils gzip
##' @importFrom R.utils gunzip
##' @importFrom downloader download
##' @export loadIntActFTP
##' @examples loadIntActFTP("./")
##' @author Vitalii Kleshchevnikov
loadIntActFTP = function(dir){
  file = paste0(dir,"intact",lastIntActRelease(), ".txt")
  file.gz = paste0(file, ".gz")
  # download MI-TAB 2.7 from IntAct ftp
  if(!file.exists(file.gz)) {
    message("... dowloading from IntAct ftp ...")
    downloader::download("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt", file)
    gzip(filename = file, destname = file.gz, remove = T, overwrite = T)
  } else {message("... loading local copy ...")}
  gunzip(filename = file.gz, destname = file, remove = F, overwrite = T)
  intact = fread(file, header = T, stringsAsFactors = F)
  unlink(file)

  result = list(data = intact, metadata = "This object contains the data from ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt and contains all molecular interaction data from the following databases: \"IntAct\", \"MINT\", \"DIP\", \"bhf-ucl\", \"MPIDB\", \"MatrixDB\", \"HPIDb\",\"I2D-IMEx\",\"InnateDB-IMEx\", \"MolCon\", \"UniProt\", \"MBInfo\"")

  class(result) = "RAW_MItab27"
  return(result)
}
