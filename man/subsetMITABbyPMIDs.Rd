% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetMITABbyPMIDs.R
\name{subsetMITABbyPMIDs}
\alias{subsetMITABbyPMIDs}
\title{Filter molecular interactions by Pubmed ID}
\usage{
subsetMITABbyPMIDs(MITABdata, PMIDs = NULL, inverse_filter = F)
}
\arguments{
\item{MITABdata}{object of any clean_MItab class (the class that contains "clean_MItab" in it's name and is initially produced by  \code{\link[PItools]{cleanMITAB}})}

\item{PMIDs}{character or character vector, return only interactions reported by these Pubmed IDs}

\item{inverse_filter}{logical, inverse filtering criteria}
}
\value{
return only interactions reported by \code{PMIDs}
}
\description{
filter molecular interaction data by Pubmed ID of publications
}
\examples{
{
# Download all human interaction data
full = fullInteractome(taxid = "9606", database = "IntActFTP",
         clean = TRUE, protein_only = TRUE)

# subset unpublished Vidal group data
Vidal = subsetMITABbyPMIDs(MITABdata = full,
               PMIDs = "unassigned1304")

# subset both published and unpublished Vidal group data
Vidal_all = subsetMITABbyPMIDs(MITABdata = full,
               PMIDs = c("25416956", "unassigned1304"))

# subset Mattias Mann 2015 paper data
Mann = subsetMITABbyPMIDs(MITABdata = full,
               PMIDs = "26496610")
}
}
\seealso{
\code{\link[PItools]{subsetMITABbyMethod}}), \code{\link[PItools]{subsetMITABbyID}})
}
\author{
Vitalii Kleshchevnikov
}
