##' 2D bin plot of the features of top N domain-protein pairs as identified by \code{\link[MItools]{permutationPval}}
##' @name PermutResult2D
##' @author Vitalii Kleshchevnikov
##' @param res \code{\link[MItools]{permutationPval}} output of (class XYZinteration_XZEmpiricalPval)
##' @param N the number of top pairs to look at
##' @param id.cols character, columns that define domain-protein pairs (domain id and protein id)
##' @param value.cols which values to plot
##' @param top_by column name by which to select top N pair
##' @return ggplot
##' @import data.table
##' @import GGally
##' @importFrom ggplot2 theme_light
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 geom_bin2d
##' @importFrom ggplot2 scale_fill_gradient
##' @importFrom ggplot2 scale_y_log10
##' @importFrom ggplot2 scale_x_log10
##' @importFrom ggplot2 annotation_logticks
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_density
##' @export PermutResult2D

PermutResult2D = function(res, N, id.cols = c("IDs_interactor_viral", "IDs_domain_human"), value.cols = c("domain_count", "IDs_interactor_viral_degree", "domain_count_per_IDs_interactor_viral", "p.value"), rank_by = "p.value", filter = NULL){

  res_temp = unique(res$data_with_pval[, unique(c(id.cols, value.cols, rank_by)), with = F])
  if(!is.null(filter)){
    res_temp = res_temp[eval(formula(paste0("~",filter))[[2]]),]
  }

  GGally::ggpairs(res_temp[order(eval(as.formula(paste0("~", rank_by))[[2]]), decreasing = F)[1:N],],
                  columns = value.cols,
                  lower = list(continuous = GGally_d2_bin_log10)#,
                  #diag = list(continuous = geom_density)
  ) +
    theme_light() +
    theme(strip.text.y = element_text(angle = 0, size = 10),
          strip.text.x = element_text(angle = 90, size = 10))
}
# function to accomodate ggplot2::geom_bin2d in GGally::ggpairs, taken from http://ggobi.github.io/ggally/#custom_functions
GGally_d2_bin_log10 <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    scale_fill_gradient(low = low, high = high) +
    scale_y_log10() + scale_x_log10() + annotation_logticks()
}

GGally_log10_density = function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_density(...) +
    scale_x_log10()
}
