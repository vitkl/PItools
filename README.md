# R Package - MItools
## Molecular interaction data tools: retrieving data from IntAct and most other databases using PSICQUIC service

Package contains high-level functions for retrieving molecular interaction data from most public databases using PSICQUIC service (http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml) and manipulating that data. Useful features include: retrieve full interactome of a given species, retrieve all interactions between two taxonomy groups (taxonomy structure-aware search), query molecular interaction databases using MIQL, subset interaction data by list of protein/gene identifiers, interaction detection method, publication (PMIDs).
    Depends on PSICQUIC, data.table.  

## Installation

```r
# Install R MItools package
install.packages("BiocManager") # for installing BioConductor dependencies
BiocManager::install("vitkl/MItools", dependencies = T)
```

## Basic use

Load all interactions between human proteins.

```r
human = fullInteractome(taxid = 9606, database = "IntActFTP",
                        format = "tab27", clean = TRUE,
                        protein_only = TRUE,
                        directory = NULL) # keep data files inside R library, 
                                          # or specify your directory 
```

Follow this [example](https://vitkl.github.io/MItools/articles/articles/Full_interactomes_interspecies.html) for more details.  

Beware of [study bias](https://vitkl.github.io/MItools/articles/articles/study_bias_doc.html) (ascertainment bias): in aggregated data better studied proteins have more interactions (not a problem for unbiased proteome-wide interaction screens).  
