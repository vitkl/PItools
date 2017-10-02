##' Plot the enrichment of domains that are likely to be mediating interactions in domain annotated by ELM as mediating linear motif mediated interactions
##' @name plotEnrichment
##' @author Vitalii Kleshchevnikov
##' @param ... objects (matrix) the output of \code{runningTestEnrichment}
##' @param list list of object (alternative way to pass the output of \code{runningTestEnrichment})
##' @param random_domains the result of comparing random domain-protein pairs (by \code{random_domains} function), can be null
##' @param domains_known_mapped domains from ELM mapped to InterProID
##' @param type what to plot: "pval"- Fisher test p-val (in out set vs in ELM), "odds_ratio" - Fisher test odds ratio, "count" of domains in our set that are in ELM
##' @param plot_type r-base plot type default is "l"
##' @param plot_name to be passed to plot(main = *)
##' @param plot_args other arguments to plot: to be passed as vector of characters, such as c("cex = 1.2", "family = \"Helvetica\"")
##' @param legend_args other arguments to legend: to be passed as vector of characters, such as c("cex = 1.2")
##' @return plot
##' @import data.table
##' @export plotEnrichment
plotEnrichment = function(..., runningTestEnrichmentlist = list(), random_domains = NULL, domains_known_mapped, type = "count", plot_type = "l", plot_name = "", plot_args = NULL, legend_args = NULL, leg_pos_x = NULL, show_known_domains = F, plot_total_domains_found = F){

  res = list(...)
  res = c(res, runningTestEnrichmentlist)
  typenum = match(type, c("pval", "odds_ratio", "count"))
  ngroups = length(res)

  if(type == "count" & show_known_domains) color = colorRampPalette(brewer.pal(7, "Dark2"))(ngroups + 1) else color = colorRampPalette(brewer.pal(7, "Dark2"))(ngroups)
  if(is.na(typenum)) stop("'type' should be one of “count”, “odds_ratio”, “pval”")

  leg_pos_y = max(sapply(res, function(x, typenum) max(as.numeric(x[typenum,])), typenum))
  if(!is.null(random_domains)) leg_pos_y = max(leg_pos_y, random_domains[typenum][[1]])
  if(is.null(leg_pos_x)){
    leg_pos_x = max(sapply(res, function(x) max(as.numeric(colnames(x))))) * 0.20
  }

  xlim_up = max(sapply(res, function(x) max(as.numeric(colnames(x)))))

  if(type == "pval") {ylim = c(0, 1); ylab = "p-value"}
  if(type == "count") {
    if(show_known_domains){
      leg_pos_y = length(domains_known_mapped) - 1
      ylim = c(0, length(domains_known_mapped) + 1)
    } else {
      ylim = c(0,leg_pos_y)
      }
    ylab = "known domain found"
    }
  if(type == "odds_ratio") {ylim = c(0,leg_pos_y); ylab = "odds ratio: in top N pairs AND in ELM / \nnot in top N pairs AND in ELM"}

  if(is.null(plot_args)){
    plot(colnames(res[[1]]),rep(0,ncol(res[[1]])),
         ylab = ylab, xlab = "top N viral protein - domain pairs selected",
         type = plot_type, ylim = ylim, lwd = 0,
         main = plot_name, xlim = c(0, xlim_up + 20))
  } else {
    additional_plot_args = paste0(plot_args, collapse = ",")
    plot_text = paste0("plot(colnames(res[[1]]),rep(0,ncol(res[[1]])),
       ylab = ylab, xlab = \"top N viral protein - domain pairs selected\",
           type = plot_type, ylim = ylim, lwd = 0,
           main = plot_name, xlim = c(0, xlim_up + 20),",additional_plot_args, ")")
    eval(parse(text = plot_text))
  }

  # plot random domains quantiles
  if(!is.null(random_domains)){
    random_legend = c("97.5% quantile", "75% quantile", "median", "25% quantile", "2.5% quantile")
    random_cols = c("#DDDDDD", "#CCCCCC", "#AAAAAA", "#CCCCCC", "#DDDDDD")
    random_line_width = c(2,4,8,4,2)
    for (i in 1:5) {
      lines(x = colnames(random_domains[typenum][[1]]), y = random_domains[typenum][[1]][random_legend[i],], col = random_cols[i], lwd = random_line_width[i], type = plot_type)
    }
  }

  for (i in 1:ngroups) {
    lines(x = colnames(res[[i]]), y = res[[i]][typenum,], col = color[i], type = plot_type, lwd = 3)
    if(type == "count" & plot_total_domains_found){
      lines(x = colnames(res[[i]]), y = res[[i]]["total_count",], col = color[i], type = plot_type, lwd = 3)
    }
  }

  if(type == "count") {
    if(show_known_domains){
    abline(h = length(domains_known_mapped), col = color[ngroups + 1])
    }
    }

  legend_names = c("statictic used in permutation test:")
  for(i in 1:ngroups){
    legend_names = c(legend_names, unique(res[[i]]["name",]))
  }

  if(type == "count" & show_known_domains) legend_names = c(legend_names, "domains known to interact with linear motifs")

  line_width = rep(3, length(color) + 1)

  if(!is.null(random_domains)){
    legend_names = c(legend_names, paste0("N random protein-domain pairs, ", random_legend))
    color = c(color, random_cols)
    line_width = c(line_width, random_line_width)
  }

  if(is.null(legend_args)){
    legend(x = leg_pos_x, y = leg_pos_y, legend_names,
           col = c("white", color), lty = 1, lwd = line_width, merge = TRUE)
  } else {
    additional_legend_args = paste0(legend_args, collapse = ",")
    legend_text = paste0("legend(x = leg_pos_x, y = leg_pos_y, legend_names,
         col = c(\"white\", color), lty = 1, lwd = line_width, merge = TRUE,",additional_legend_args, ")")
    eval(parse(text = legend_text))
  }

}
