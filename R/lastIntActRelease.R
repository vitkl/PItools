##' last IntAct (all IMEx databases) release
##' @name lastIntActRelease
##' @return character string specifying the year, month, day, and the hour when the IntAct data file was last modified
##' @import data.table
##' @export lastIntActRelease
##'
lastIntActRelease = function(){
  last_release = fread("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/")
  last_release = last_release[V9 == "all.zip", paste0(V6,V7, "_", gsub(":",".", V8))]
  return(paste0(year(Sys.Date()), last_release))
}
