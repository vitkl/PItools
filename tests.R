library(roxygen2)
library(devtools)
document()
install()

library(MItools)
library(GenomicRanges)
library(R.utils)
data = fread("../viral_project/processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)
file_BioPlex3 = "../viral_project/processed_data_files/human_net_w_domains"
gunzip(paste0(file_BioPlex3,".gz"), remove = F)
data = fread(file_BioPlex3, sep = "\t", stringsAsFactors = F)
data = data[IDs_interactor_human_B != "UNKNOWN"]
unlink(file_BioPlex3)
res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                      associations2test = IDs_interactor_viral ~ IDs_domain_human,
                      node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                       IDs_domain_human ~ domain_count,
                                       IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                      data = data,
                      statistic = IDs_interactor_viral + IDs_domain_human ~ .N,
                      select_nodes = IDs_domain_human ~ domain_count >= 1,
                      N = 100, clustermq_jobs = 20, split_comp_inner_N = 2,
                      cores = NULL, seed = 2, clustermq = T, clustermq_mem = 2000)
plot(res)
res

microbenchmark::microbenchmark({res <- permutationPval(interactions2permute = IDs_interactor_human_A ~ IDs_interactor_human_B, # first set of interacting pairs (XY) that are to be permuted
                                                       associations2test = IDs_interactor_human_A ~ IDs_domain_human_B, # set of interacting pairs to be tested (XZ), YZ interactions are assumed
                                                       node_attr = list(IDs_interactor_human_A ~ IDs_interactor_human_A_degree,
                                                                        IDs_domain_human_B ~ domain_count),
                                                       data = data,
                                                       statistic = IDs_interactor_human_A + IDs_domain_human_B ~ .N,
                                                       select_nodes = IDs_domain_human_B ~ domain_count >= 1,
                                                       N = 180,
                                                       cores = NULL, seed = 2, also_permuteYZ = F,
                                                       clustermq = T, clustermq_mem = 20000,
                                                       split_comp_inner_N = 3, clustermq_jobs = 30,
                                                       clustermq_log_worker = F)}, times = 1)
plot(res, IDs_interactor_human_A + IDs_domain_human_B ~ log10(not_missing))
res$data_with_pval[!is.na(IDs_domain_human_B)]

{permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                 associations2test = IDs_interactor_viral ~ IDs_domain_human,
                 node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                  IDs_domain_human ~ domain_count,
                                  IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                 data = data,
                 statistic = IDs_interactor_viral + IDs_domain_human ~ .N,
                 select_nodes = IDs_domain_human ~ domain_count >= 1,
                 N = 1000,
                 cores = 3, seed = 2, clustermq = F)}


set.seed(1)
random = randomInteractome(n_prot = 200, degree_dist = NULL, taxid = "9606", database = "imex", protein_only = TRUE)

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
                       cores = NULL, seed = NULL)

microbenchmark::microbenchmark({res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                                                      associations2test = IDs_interactor_viral ~ IDs_domain_human,
                                                      node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                                                       IDs_domain_human ~ domain_count,
                                                                       IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                                                      data = data,
                                                      statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                                                      select_nodes = IDs_domain_human ~ domain_count >= 1,
                                                      N = 10,
                                                      cores = NULL, seed = 1)}, times = 10)

profvis::profvis({resEnv = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                                           associations2test = IDs_interactor_viral ~ IDs_domain_human,
                                           node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                                            IDs_domain_human ~ domain_count,
                                                            IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                                           data = data,
                                           statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                                           select_nodes = IDs_domain_human ~ domain_count >= 1,
                                           N = 10,
                                           cores = NULL, seed = 1)})

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
                                                            cores = NULL, seed = 1)}, times = 10)

qplot(x = resFISHER$data_with_pval[p.value < 0.5, IDs_interactor_viral_degree], y = resFISHER$data_with_pval[p.value < 0.5, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()
qplot(x = res$data_with_pval[p.value < 0.5, IDs_interactor_viral_degree], y = res$data_with_pval[p.value < 0.5, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()

qplot(x = resFISHER$data_with_pval[p.value < 0.01, IDs_interactor_viral_degree], y = resFISHER$data_with_pval[p.value < 0.01, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()
qplot(x = res$data_with_pval[p.value < 0.01, IDs_interactor_viral_degree], y = res$data_with_pval[p.value < 0.01, domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()

# Unit: seconds (without inner and outer replicate)
#expr
#{     resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~          IDs_interactor_human, associations2test = IDs_interactor_viral ~          IDs_domain_human, node_attr = list(IDs_interactor_viral ~          IDs_interactor_viral_degree, IDs_domain_human ~ domain_count +          N_prot_w_interactors, IDs_interactor_viral + IDs_domain_human ~          domain_count_per_IDs_interactor_viral), data = data,          statistic = IDs_interactor_viral + IDs_domain_human ~              fisher.test(matrix(c(unique(domain_count), unique(N_prot_w_interactors) -                  unique(domain_count), unique(domain_count_per_IDs_interactor_viral),                  unique(IDs_interactor_viral_degree) - unique(domain_count_per_IDs_interactor_viral)),                  2, 2), alternative = "greater", conf.int = F)$p.value,          select_nodes = IDs_domain_human ~ domain_count >= 1,          N = 100, cores = NULL, seed = 1) }
#min       lq    mean   median       uq      max neval
#34.89687 35.35519 35.3463 35.39729 35.42503 35.49184    10

# Unit: seconds (with inner and outer replicate)
# expr
# {     resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~          IDs_interactor_human, associations2test = IDs_interactor_viral ~          IDs_domain_human, node_attr = list(IDs_interactor_viral ~          IDs_interactor_viral_degree, IDs_domain_human ~ domain_count +          N_prot_w_interactors, IDs_interactor_viral + IDs_domain_human ~          domain_count_per_IDs_interactor_viral), data = data,          statistic = IDs_interactor_viral + IDs_domain_human ~              fisher.test(matrix(c(unique(domain_count), unique(N_prot_w_interactors) -                  unique(domain_count), unique(domain_count_per_IDs_interactor_viral),                  unique(IDs_interactor_viral_degree) - unique(domain_count_per_IDs_interactor_viral)),                  2, 2), alternative = "greater", conf.int = F)$p.value,          select_nodes = IDs_domain_human ~ domain_count >= 1,          N = 100, cores = NULL, seed = 1) }
# min       lq    mean   median       uq      max neval
# 38.02992 38.33805 38.7546 38.73306 38.86937 39.60011    10


library(MItools)
library(rtracklayer)
library(ggplot2)
data = fread("../viral_project/processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)
time = proc.time()
res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                      associations2test = IDs_interactor_viral ~ IDs_domain_human,
                      node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                       IDs_domain_human ~ domain_count,
                                       IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                      data = data,
                      statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                      select_nodes = IDs_domain_human ~ domain_count >= 1,
                      N = 10000,
                      cores = NULL, seed = 2)
proc.time() - time
time = proc.time()
resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                            associations2test = IDs_interactor_viral ~ IDs_domain_human,
                            node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                             IDs_domain_human ~ domain_count + N_prot_w_interactors,
                                             IDs_interactor_viral + IDs_domain_human ~ domain_count_per_IDs_interactor_viral),
                            data = data,
                            statistic = IDs_interactor_viral + IDs_domain_human ~ fisher.test(matrix(c(unique(domain_count), unique(N_prot_w_interactors) - unique(domain_count), unique(domain_count_per_IDs_interactor_viral), unique(IDs_interactor_viral_degree) - unique(domain_count_per_IDs_interactor_viral)),2,2), alternative = "greater", conf.int = F)$p.value,
                            select_nodes = IDs_domain_human ~ domain_count >= 1,
                            N = 10000,
                            cores = NULL, seed = 1)
resFISHER$data_with_pval[, p.value := 1 - p.value]
proc.time() - time

resFISHER$IDs_interactor_viral_degreeVSdomain_count = qplot(x = resFISHER$data_with_pval[order(p.value)[1:250], IDs_interactor_viral_degree], y = resFISHER$data_with_pval[order(p.value)[1:250], domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()
res$IDs_interactor_viral_degreeVSdomain_count = qplot(x = res$data_with_pval[order(p.value)[1:250], IDs_interactor_viral_degree], y = res$data_with_pval[order(p.value)[1:250], domain_count], geom = "bin2d") + scale_x_log10() + scale_y_log10()

interactiondomains = fread("http://elm.eu.org/interactiondomains.tsv")
interactiondomains[, pfam_id := `Interaction Domain Id`]

domains_known = interactiondomains[, unique(pfam_id)]

"../viral_project/processed_data_files/InterProScan_domains_nonredundant.gff3"
"../viral_project/processed_data_files/all_human_viral_protein_domains.gff3.gz"
InterProScan_domains_nonred = import(con = "../viral_project/processed_data_files/InterProScan_domains_nonredundant.gff3", format = "gff3")
domains_mapping = unique(data.table(any_id = as.character(InterProScan_domains_nonred$Name), interpro_id = as.character(InterProScan_domains_nonred$Dbxref)))
domains_known_mapped = unique(domains_mapping[any_id %in% domains_known, interpro_id])
domains_not_mapped = unique(domains_known[!domains_known %in% domains_mapping$any_id])

test_enrichment = function(N, res, domains_known_mapped){
  res$data_pval = unique(res$data_with_pval[,.(IDs_interactor_viral, IDs_domain_human, p.value, domain_type, domain_count, IDs_interactor_viral_degree)])
  res$data_pval[, pval_fdr := p.adjust(p.value, method = "fdr")]
  hist(res$data_pval[, pval_fdr], breaks = seq(0,1,0.01))

  domains_found = res$data_pval[order(p.value)[1:N], unique(IDs_domain_human)]

  alldomains = res$data_pval[, unique(IDs_domain_human)]
  known = factor(alldomains %in% domains_known_mapped, levels = c("TRUE", "FALSE"))
  found = factor(alldomains %in% domains_found, levels = c("TRUE", "FALSE"))
  table_res = table(known, found)

  test = fisher.test(table(known, found), alternative = "greater", conf.int = T)

  return(c(test$p.value, test$estimate, table_res["TRUE", "TRUE"]))
}

enrichment = sapply(seq(25, 500, 25), test_enrichment, res, domains_known_mapped)
colnames(enrichment) = seq(25, 500, 25)
enrichmentFISHER = sapply(seq(25, 500, 25), test_enrichment, resFISHER, domains_known_mapped)
colnames(enrichmentFISHER) = seq(25, 500, 25)

plot(colnames(enrichment), enrichment[2,], ylab = "Fisher test odds ratio", xlab = "top N viral protein - domain pairs selected", col = "red", type = "l", ylim = c(0,18))
lines(x = colnames(enrichment), y = enrichmentFISHER[2,], col = "blue", type = "l")
legend(x = 80, y = 17.5, c("statictic used in permutation test:","domain frequency among interactors of a viral protein", "Fisher test pval: domain overrepresentation over the background"), col = c("white","red", "blue"), lty = 1 ,merge = TRUE)

plot(colnames(enrichment), enrichment[3,], ylab = "known domain found", xlab = "top N viral protein - domain pairs selected", col = "red", type = "l", ylim = c(0,length(domains_known_mapped)+1))
lines(x = colnames(enrichment), y = enrichmentFISHER[3,], col = "blue", type = "l")
abline(h = length(domains_known_mapped), col = "green")
legend(x = 80, y = 50, c("statictic used in permutation test:","domain frequency among interactors of a viral protein", "Fisher test pval: domain overrepresentation over the background", "domains known to interact with linear motifs"), col = c("white","red", "blue", "green"), lty = 1 , merge = TRUE)

plot(colnames(enrichment), enrichment[1,], ylab = "Fisher test pvalue", xlab = "top N viral protein - domain pairs selected", col = "red", type = "l", ylim = c(0,0.004))
lines(x = colnames(enrichment), y = enrichmentFISHER[1,], col = "blue", type = "l")
legend(x = 80, y = 0.0041, c("statictic used in permutation test:","domain frequency among interactors of a viral protein", "Fisher test pval: domain overrepresentation over the background"), col = c("white","red", "blue"), lty = 1 ,merge = TRUE)

big_jobs = sapply(list.files(), function(file){length(readLines(file))}) == 1
error_paths = paste0("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/qslimfinder.Full_IntAct4.FALSE/log_dir/error/", gsub("\\.sh","", names(big_jobs)[big_jobs]))
sapply(error_paths, function(error_path) {
  system(paste0("cat ", error_path," | grep Terminated"), intern=T)
})

log_paths = paste0("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/qslimfinder.Full_IntAct4.FALSE/log_dir/log/", gsub("\\.sh","", names(big_jobs)[big_jobs]))
job_status = sapply(log_paths, function(log_path) {
  system(paste0("cat ", log_path," | grep TERM_MEMLIMIT"), intern=T)
})

# how many terminated because of memory
sum(job_status == "TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.") # 9
sum(sapply(job_status, function(element) {
  length(element) == 0
})) # 65

all.jobs = sapply(list.files(), function(file){length(readLines(file))})
error_paths = paste0("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/qslimfinder.Full_IntAct4.FALSE/log_dir/error/", gsub("\\.sh","", names(all.jobs)))
errors = sapply(error_paths, function(error_path) {
  system(paste0("cat ", error_path," | grep Terminated"), intern=T)
})

log_paths = paste0("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/qslimfinder.Full_IntAct4.FALSE/log_dir/log/", gsub("\\.sh","", names(all.jobs)))
successes = sapply(log_paths, function(log_path) {
  system(paste0("cat ", log_path," | grep \"Successfully completed\""), intern=T)
})

## TEST for benchmarking
occurence_file = "../viral_project/qslimfinder.Full_IntAct3.FALSE/result/occurence.txt"
main_file = "../viral_project/qslimfinder.Full_IntAct3.FALSE/result/main_result.txt"
domain_res_file = "../viral_project/processed_data_files/domain_res_count_20171019.RData"
motif_setup = "../viral_project/processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Full_IntAct3.FALSE_clust201802.RData"
neg_set = c("all_instances", "all_proteins", "random")[1]
domain_results_obj = "res_count"
motif_input_obj = "forSLIMFinder_Ready"
motif_setup_obj2 = NULL
occurence_filt = NULL
one_from_cloud = T
dbfile_main = "../viral_project/data_files/instances_all.gff"
dburl_main = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic="
dbfile_query = "../viral_project/data_files/instances_query.gff"
dburl_query = "http://elm.eu.org/instances.gff?q=all&taxon=irus&instance_logic="
query_res_query_only = T
motif_types = c("DOC", "MOD", "LIG")
all_res_excl_query = T
merge_motif_variants = F
seed = 21
N = 100
replace = T
within1sequence = T
query_predictor_col = "Sig"
all_predictor_col = "Sig"
normalise = T
minoverlap = 2
minoverlap_redundant = 5
merge_domain_data = T
merge_by_occurence_mcols = c("query", "interacts_with")
merge_by_domain_res_cols = c("IDs_interactor_viral", "IDs_interactor_human", "IDs_domain_human", "Taxid_interactor_human","Taxid_interactor_viral")
merge_by_non_query_domain_res_cols = c("IDs_interactor_human_A", "IDs_interactor_human_B", "IDs_domain_human_B", "Taxid_interactor_human_A","Taxid_interactor_human_B")
count_ranges_by = list(by = "IDs_domain_human", name = "motif_occ_per_domain",
                       normalise_by = "domain_count", normalised_name = "normalised_motif_occ_per_domain")
filter_by_domain_data = "p.value < 0.5"
motif_pval_cutoff = 1
select_predictor_per_range = max
non_query_domain_res_file = "../viral_project/processed_data_files/predict_domain_human_clust20180819.RData"
non_query_domain_results_obj = "res_count_all" # NULL res_count_all
non_query_domains_N = 0
non_query_set_only = F
query_domains_only = F
min_non_query_domain_support = 3
select_top_domain = F

# creating a website
# Install release version from CRAN
install.packages("pkgdown")
pkgdown::build_site()


## benchmarking ways to combine p-values for domains

prior = mcol[p.value<1,mean(1-p.value)]
mcol[, combined_p.value := updatePvalue(prior, p.value), by = .(query, Pattern, interacts_with, domain)]

dom_elm = merge(mcol, elm_interactions,
                by.x = c("query", "interacts_with", "domain"),
                by.y = c("interactorElm", "interactorDomain", "InterProID"),
                all.x = T, all.y = F)
dom_elm[!is.na(Elm), Correct_SLiM_binding := "yes"]
dom_elm[is.na(Elm), Correct_SLiM_binding := "no"]

# filter and select unique to have one row per protein_with_motif - domain pair
dom3 = unique(dom_elm[,.(query, domain, p.value = combined_p.value, Correct_SLiM_binding)])
dom3[, both_SLIM_interaction := uniqueN(Correct_SLiM_binding) >= 2,
     by = .(query)]
dom3 = dom3[both_SLIM_interaction == T]
# calculate N domain protein pairs
dom3[, N_pairs := uniqueN(paste0(query, domain)), by = .(Correct_SLiM_binding)]
dom3[, network := paste0(unique(N_pairs[Correct_SLiM_binding == "yes"]),
                         " correct domains ",
                         "/ ", unique(N_pairs[Correct_SLiM_binding == "no"]),
                         " total")]
ggplot(dom3, aes(p.value, color = Correct_SLiM_binding, fill = Correct_SLiM_binding)) +
  geom_histogram(alpha = 0.3, aes(y = ..density..), bins = 50, position = "identity") +
  theme_bw() +
  facet_grid(network~., scales = "free_y")+
  theme(legend.position = "none",
        strip.text.y = element_text(angle = 0),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        text = element_text(size = 14))
### benchmark end

prior = mcol[p.value<1,mean(1-p.value)]
prior
updatePvalue = function(p.value1, p.value2, N = 1, incr = 1e-10){
  N = max(length(p.value1), length(p.value2))
  # add/substract small value to p-value to avoid 0 likelihood
  p.value2[p.value2 == 1] = 1 - incr
  p.value2[p.value2 == 0] = incr
  # calculate likelihood
  likelihood_of_2_1 = dbeta(p.value2, p.value1 * N, (1 - p.value1) * N)
  likelihood_of_2_not_1 = dbeta(p.value2, (1 - p.value1) * N, p.value1 * N)

  p.value = sum(likelihood_of_2_1 * p.value1) /
    (sum(likelihood_of_2_1 * p.value1) +
       sum(likelihood_of_2_not_1 * (1 - p.value1)))
  p.value
}

probs = c(prior, 0.999, 0, 0.979)
b = updatePvalue(probs[1], probs[2:4], N = 1)
b
up2 = updatePvalue(up1$p.value, probs[3], N = up1$N)
up3 = updatePvalue(up2$p.value, probs[4], N = up2$N)
up3

combn(c(0.999, 0, 0.979),2)

fit = beta.mle(c(0.999, 0.998, 0.979), tol = 1e-09)
hist(rbeta(10000, 0.4999 + 100, 0.4999 + 100))
mean(fit$param[1] / (fit$param[1] + fit$param[2]))

### draft end
