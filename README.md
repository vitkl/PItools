# R Package - PItools
## Protein interaction data tools: retrieving data from IntAct and most other databases using PSICQUIC service

Package contains high-level functions for retrieving molecular interaction data from most public databases using PSICQUIC service (http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml) and manipulating that data. Useful features include: retrieve full interactome of a given species, retrieve all interactions between two taxonomy groups (taxonomy structure-aware search), query molecular interaction databases using MIQL, subset interaction data by list of protein/gene identifiers, interaction detection method, publication (PMIDs).
    Depends on PSICQUIC, data.table.  

## Installation

```r
# Install R PItools package
install.packages("BiocManager") # for installing BioConductor dependencies
BiocManager::install("Biostrings", "remotes") # dependency for PSICQUIC & for installing from github
BiocManager::install("vitkl/PItools", dependencies = T)
```

## Basic use

Load all interactions between human proteins in IntAct database.

Follow this [example](https://vitkl.github.io/PItools/articles/articles/Full_interactomes_interspecies.html) for more details.  

```r
# load package
library(PItools)

# load all interactions from IntAct FTP storage (https://www.ebi.ac.uk/intact/downloads)
human = fullInteractome(taxid = 9606, database = "IntActFTP",
                        format = "tab27",
                        clean = TRUE, # parse into usable format (takes 5-10 minutes)
                        protein_only = TRUE, # filter protein interactions
                        directory = NULL) # keep data files inside R library, 
                                          # or specify your directory 
```
```r
# protein interactions lack directionality & are detected in multiple studies
# so find the number of unique interactions
uniqueNinteractions(human)
uniqueNinteractors(human)
```
```r
# subset by a list of proteins
mediator_complex_proteins = fread("https://www.uniprot.org/uniprot/?query=GO:0016592%20AND%20taxonomy:9606&format=tab&columns=id")
mediator_complex = subsetMITABbyID(human, ID_seed = mediator_complex_proteins$Entry,
                                   within_seed = T, only_seed2nonseed = F)
uniqueNinteractions(mediator_complex)
uniqueNinteractors(mediator_complex)
```

Load interactions for one protein identified with a specific method.

```r
# Query for interactions of bacterial RNA polymerase sigma factor SigA 
# identified using two-hybrid methods inall imex databases 
# (result identical to IntActFTP but only requested interactions are downloaded)
inter = queryPSICQUICrlib(query = "id:P74565 AND detmethod:\"MI:0018\"",
                  format = "tab25",
                  database = "imex",
                  directory = "./data_files/")
```
```r
# The data returned by queryPSICQUICrlib constains auxillary information 
# that is not necessary for most analysis. Let's clean the data.
inter = cleanMITAB(inter)
inter$data[1:5]
```

Beware of [study bias](https://vitkl.github.io/PItools/articles/articles/study_bias_doc.html) (ascertainment bias): in aggregated data better studied proteins have more interactions (not a problem for unbiased proteome-wide interaction screens).  
