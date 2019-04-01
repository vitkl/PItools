##' \code{foldEnrichment} calculates fold enrichment of a domain (or other feature) among interacting partners of a protein
##' @name foldEnrichment
##' @param net data.table that specifies the PPI network:
##' 1. first column specifies a protein ("IDs_interactor_viral")
##' 2. second column specifies it's interacting partner which annotations are of interest ("IDs_interactor_human")
##' 3. third column specifies the degree of a protein in the first column ("IDs_interactor_viral_degree")
##' @param protein_annot data.table that specifies annotations for a protein in the second column of \code{net}
##' 1. first column specifies a protein ("IDs_interactor_human")
##' 2. second column specifies features of that protein ("IDs_domain_human")
##' 3. third column specifies background frequency of those features ("domain_frequency)
##' @param frequency fold enrichment or frequency in a set (if TRUE - frequency)
##' @return data.table containing fold enrichment or frequency for each domain - protein pair. Columns: IDs_interactor_viral, IDs_interactor_human, IDs_domain_human, domain_frequency_per_set or fold_enrichment
##' @author Vitalii Kleshchevnikov
##' @import data.table

foldEnrichment = function(net, protein_annot, frequency = T){
  if(ncol(net) != 3 | mean(c("IDs_interactor_viral", "IDs_interactor_human", "IDs_interactor_viral_degree") %in% colnames(net)) != 1) stop("net contains more or less columns than required or wrong colnames")
  if(ncol(protein_annot) != 3 | mean(c("IDs_interactor_human", "IDs_domain_human", "domain_frequency") %in% colnames(protein_annot)) != 1) stop("protein_annot contains more or less columns than required or wrong colnames")

  # add domain annotation to the network, keep proteins without domains (delete domains not corresponding to the proteins in the network)
  merged_net = unique(protein_annot[net, on = "IDs_interactor_human", allow.cartesian = T])
  merged_net[IDs_domain_human == "", IDs_domain_human := NA]

  # count human proteins with specific domain per viral protein
  # (how many proteins the domain is located in) per viral protein (ID) and human domain (ID)
  merged_net[, domain_count_per_IDs_interactor_viral := .N, by = .(IDs_interactor_viral, IDs_domain_human)]
  # 0 domain viral protein
  merged_net[is.na(IDs_domain_human), domain_count_per_IDs_interactor_viral := 0]

  # domain frequency but per viral protein
  merged_net[, domain_frequency_per_set := domain_count_per_IDs_interactor_viral / IDs_interactor_viral_degree]
  # 0 domain viral protein
  merged_net[is.na(IDs_domain_human), domain_frequency_per_set := 0]
  if(frequency) merged_net = merged_net[,.(IDs_interactor_viral, IDs_interactor_human, IDs_domain_human, domain_frequency_per_set)]
  if(!frequency){
    # fold enrichment
    merged_net[, fold_enrichment := domain_frequency_per_set / domain_frequency]
    # 0 domain viral protein
    merged_net[is.na(IDs_domain_human), fold_enrichment := 0]

    merged_net = merged_net[,.(IDs_interactor_viral, IDs_interactor_human, IDs_domain_human, fold_enrichment)]
  }
  return(merged_net)
}
