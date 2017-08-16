library(roxygen2)
library(devtools)
document()
install()

library(MItools)
data = fread("../viral_project/processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)
res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                associations2test = IDs_interactor_viral ~ IDs_domain_human,
                node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                 IDs_domain_human ~ domain_count,
                                 IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                data = data,
                statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                select_nodes = IDs_domain_human ~ domain_count >= 1,
                N = 100,
                cores = NULL, seed = 2, include_missing_Z_as_zero = T)

res_g = res
all.equal(res_g, res)

res2 = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                      associations2test = IDs_interactor_viral ~ IDs_domain_human,
                      node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                       IDs_domain_human ~ domain_count,
                                       IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                      data = data,
                      statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                      select_nodes = IDs_domain_human ~ domain_count > 16,
                      N = 10,
                      cores = NULL, seed = NULL, include_missing_Z_as_zero = T)

microbenchmark::microbenchmark({res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                                                      associations2test = IDs_interactor_viral ~ IDs_domain_human,
                                                      node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                                                       IDs_domain_human ~ domain_count,
                                                                       IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                                                      data = data,
                                                      statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                                                      select_nodes = IDs_domain_human ~ domain_count >= 1,
                                                      N = 10,
                                                      cores = NULL, seed = 1, include_missing_Z_as_zero = T)}, times = 10)

profvis::profvis({resEnv = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                                        associations2test = IDs_interactor_viral ~ IDs_domain_human,
                                        node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                                         IDs_domain_human ~ domain_count,
                                                         IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                                        data = data,
                                        statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                                        select_nodes = IDs_domain_human ~ domain_count >= 1,
                                        N = 10,
                                        cores = NULL, seed = 1, include_missing_Z_as_zero = T)})

# Fisher test
microbenchmark::microbenchmark({resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                                                      associations2test = IDs_interactor_viral ~ IDs_domain_human,
                                                      node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                                                       IDs_domain_human ~ domain_count + N_prot_w_interactors,
                                                                       IDs_interactor_viral + IDs_domain_human ~ domain_count_per_IDs_interactor_viral),
                                                      data = data,
                                                      statistic = IDs_interactor_viral + IDs_domain_human ~ fisher.test(matrix(c(unique(domain_count), unique(N_prot_w_interactors) - unique(domain_count), unique(domain_count_per_IDs_interactor_viral), unique(IDs_interactor_viral_degree) - unique(domain_count_per_IDs_interactor_viral)),2,2), alternative = "greater", conf.int = F)$p.value,
                                                      select_nodes = IDs_domain_human ~ domain_count >= 1,
                                                      N = 100,
                                                      cores = NULL, seed = 1, include_missing_Z_as_zero = T)}, times = 10)

qplot(x = resFISHER$data_with_pval[p.value < 0.5, IDs_interactor_viral_degree], y = resFISHER$data_with_pval[p.value < 0.5, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()
qplot(x = res$data_with_pval[p.value < 0.5, IDs_interactor_viral_degree], y = res$data_with_pval[p.value < 0.5, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()

qplot(x = resFISHER$data_with_pval[p.value < 0.01, IDs_interactor_viral_degree], y = resFISHER$data_with_pval[p.value < 0.01, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()
qplot(x = res$data_with_pval[p.value < 0.01, IDs_interactor_viral_degree], y = res$data_with_pval[p.value < 0.01, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()

# Unit: seconds (without inner and outer replicate)
#expr
#{     resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~          IDs_interactor_human, associations2test = IDs_interactor_viral ~          IDs_domain_human, node_attr = list(IDs_interactor_viral ~          IDs_interactor_viral_degree, IDs_domain_human ~ domain_count +          N_prot_w_interactors, IDs_interactor_viral + IDs_domain_human ~          domain_count_per_IDs_interactor_viral), data = data,          statistic = IDs_interactor_viral + IDs_domain_human ~              fisher.test(matrix(c(unique(domain_count), unique(N_prot_w_interactors) -                  unique(domain_count), unique(domain_count_per_IDs_interactor_viral),                  unique(IDs_interactor_viral_degree) - unique(domain_count_per_IDs_interactor_viral)),                  2, 2), alternative = "greater", conf.int = F)$p.value,          select_nodes = IDs_domain_human ~ domain_count >= 1,          N = 100, cores = NULL, seed = 1, include_missing_Z_as_zero = T) }
#min       lq    mean   median       uq      max neval
#34.89687 35.35519 35.3463 35.39729 35.42503 35.49184    10

# Unit: seconds (with inner and outer replicate)
# expr
# {     resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~          IDs_interactor_human, associations2test = IDs_interactor_viral ~          IDs_domain_human, node_attr = list(IDs_interactor_viral ~          IDs_interactor_viral_degree, IDs_domain_human ~ domain_count +          N_prot_w_interactors, IDs_interactor_viral + IDs_domain_human ~          domain_count_per_IDs_interactor_viral), data = data,          statistic = IDs_interactor_viral + IDs_domain_human ~              fisher.test(matrix(c(unique(domain_count), unique(N_prot_w_interactors) -                  unique(domain_count), unique(domain_count_per_IDs_interactor_viral),                  unique(IDs_interactor_viral_degree) - unique(domain_count_per_IDs_interactor_viral)),                  2, 2), alternative = "greater", conf.int = F)$p.value,          select_nodes = IDs_domain_human ~ domain_count >= 1,          N = 100, cores = NULL, seed = 1, include_missing_Z_as_zero = T) }
# min       lq    mean   median       uq      max neval
# 38.02992 38.33805 38.7546 38.73306 38.86937 39.60011    10


interactiondomains = fread("http://elm.eu.org/interactiondomains.tsv")
interactiondomains[, pfam_id := `Interaction Domain Id`]

domains_known = interactiondomains[, unique(pfam_id)]

"../viral_project/processed_data_files/InterProScan_domains_nonredundant.gff3"
"../viral_project/processed_data_files/all_human_viral_protein_domains.gff3.gz"
InterProScan_domains_nonred = import(con = "../viral_project/processed_data_files/InterProScan_domains_nonredundant.gff3", format = "gff3")
domains_mapping = unique(data.table(any_id = as.character(InterProScan_domains_nonred$Name), interpro_id = as.character(InterProScan_domains_nonred$Dbxref)))
domains_known_mapped = unique(domains_mapping[any_id %in% domains_known, interpro_id])
domains_not_mapped = unique(domains_known[!domains_known %in% domains_mapping$any_id])

res$data_with_pval[, pval_fdr := p.adjust(p.value, method = "fdr")]
hist(res$data_with_pval[, pval_fdr], breaks = seq(0,1,0.01))

domains_found = res$data_with_pval[pval_fdr < 0.05, unique(IDs_domain_human)]

alldomains = res$data_with_pval[, unique(IDs_domain_human)]
known = factor(alldomains %in% domains_known_mapped, levels = c("TRUE", "FALSE"))
found = factor(alldomains %in% domains_found, levels = c("TRUE", "FALSE"))
table(known, found)
fisher.test(table(set, found), alternative = "greater", conf.int = T)


