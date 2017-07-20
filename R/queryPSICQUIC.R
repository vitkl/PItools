##' Use PSICQUIC to query for molecular interactions (using MIQL, specifying database and output format)
##' @name queryPSICQUIC
##' @author Vitalii Kleshchevnikov
##' @description query_PSICQUIC function queries a specified molecular interaction database using MIQL query language, retrieves the data in a specified format and saves it to a file
##' @details See \url{https://psicquic.github.io/MiqlReference.html} or \url{https://psicquic.github.io/MiqlReference27.html} for the description of MIQL query language.
##' @details Unlike rawQuery function from PSICQUIC package this function allows to use all IMEx databases by passing "imex" to database argument and also separates the choice of output format from the query. Output format options: \url{https://psicquic.github.io/formats.html}
##' @details This function also splits large query result into chunks (2500 interactions) to limit the load on the PSICQUIC servers.
##' @details List of data provider names accepted by the function: "imex", "APID Interactomes", "BioGrid", "bhf-ucl", "ChEMBL", "HPIDb", "InnateDB", "InnateDB-All", "IntAct", "mentha", "MPIDB", "MatrixDB", "MINT", "Reactome", "Reactome-FIs", "EBI-GOA-miRNA", "I2D", "I2D-IMEx", "InnateDB-IMEx", "MolCon", "UniProt", "MBInfo", "VirHostNet", "BAR", "EBI-GOA-nonIntAct", "ZINC"
##' @param query search query in MIQL query language, as you would type in PSICQUIC View client: \url{http://www.ebi.ac.uk/Tools/webservices/psicquic/view/home.xhtml}
##' @param format output format (most widely used tabular format is tab25, use tab27 for more data columns)
##' @param database PSICQUIC service (full list: \url{http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS}), use "imex" shorthand to choose all IMEx databases
##' @param file file name and directory to save the result ("/dir/filename")
##' @return saves query result to a file, returns the query settings, how many interactions per database retrieved and which databases are inactive
##' @import data.table
##' @import PSICQUIC
##' @export queryPSICQUIC
##' @examples queryPSICQUIC(query = "id:P74565 AND detmethod:\"MI:0018\"",
##'                format = "tab27",
##'                database = "imex",
##'                file = "P74565_2H_interactions_imex_tab27.tsv")
##'
##' @examples queryPSICQUIC(query = "id:P74565 AND detmethod:\"MI:0018\"",
##'                format = "tab25",
##'                database = "mentha",
##'                file = "P74565_2H_interactions_mentha_tab25.tsv")
##'
##' @examples queryPSICQUIC(query = "id:P74565",
##'                format = "tab25",
##'                database = "mentha",
##'                file = "P74565_2H_interactions_mentha_tab25.tsv")
##'
##' @examples queryPSICQUIC(query = "id:156",
##'                format = "tab25",
##'                database = "BioGrid",
##'                file = "entrezgene156_interactions_BioGrid_tab25.tsv")
##'
##' @author Vitalii Kleshchevnikov
queryPSICQUIC = function(query, format, database, file){
  ## Load PSICQUIC functionality
  psicquic <- PSICQUIC()
  providers <- providers(psicquic)
  if(database == "imex"){
    ## Choose IMEx databases
    databases <- c("IntAct", "MINT", "bhf-ucl", "MPIDB", "MatrixDB",
                   "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo")
  }
  if(database != "imex"){
    databases = database
  }

  query = gsub("\"","%22", query)

  SPECIES_ID_interactome = data.table()
  N_SPECIES_ID_interactome = character(length = length(databases))
  NO_DATABASE = character(length = length(databases))
  for(indices in 1:length(databases)) {
    if(databases[indices] %in% providers) {
      ## Query for the number of interactions
      PSICQUIC_query1 = paste0(query, "?format=count")
      N_interactions <- unlist(rawQuery(psicquic, databases[indices], PSICQUIC_query1))
      ## if there are any interactions - Query for the interactions by 2500 at a time
      if(N_interactions > 0){
        N_start = 1
        N_nrows = 2500
        for(n_starts in seq(from = N_start, to = N_interactions, by = N_nrows)){
          PSICQUIC_query2 = paste0(query, "?format=",format,"&firstResult=", n_starts,
                                   "&maxResults=", N_nrows)
          SPECIES_ID_interactome_d <- as.data.table(rawQuery(psicquic, databases[indices], PSICQUIC_query2))
          SPECIES_ID_interactome <- rbind(SPECIES_ID_interactome, SPECIES_ID_interactome_d)
        }
      }
      # record how many interactions per database
      N_SPECIES_ID_interactome[indices] = N_interactions
    }
    ## record if the database is not active
    else {
      NO_DATABASE[indices] = databases[indices]
    }
  }

  # add informative header for MI-TAB 2.7
  if(format == "tab27") colnames(SPECIES_ID_interactome) = unlist(strsplit(readLines("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt",  n = 1), "\t"))

  # write results to table
  if(nrow(SPECIES_ID_interactome) > 0) fwrite(SPECIES_ID_interactome, file, sep = "\t") else message("no interactions matching your query")

  # returns the query settings, how many interactions per database retrieved and which databases are inactive
  res_summary = data.table(query = query,
                           file = file,
                           format = format,
                           all.databases = databases,
                           n.interactions.in.database = N_SPECIES_ID_interactome,
                           database.not.active = NO_DATABASE)
  return(res_summary)
}
