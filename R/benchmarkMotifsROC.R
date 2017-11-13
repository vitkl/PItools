##' Read QSLIMFinder (SLIMFinder) occurence output to Genomic Ranges object
##' @rdname benchmarkMotifsROC
##' @name benchmarkMotifsROC
##' @author Vitalii Kleshchevnikov
##' @param res object class \code{benchmarkMotifsResult}, the output of \code{\link{benchmarkMotifs}}
##' @param type benchmark (test if known motifs were observed) in query proteins or in all proteins the sequences of which were used in QSLIMFinder search
##' @param ... args passed to \code{\link[ROCR]{plot}}
##' @param text_args character 1L, arguments to \code{\link[base]{text}} and \code{\link[base]{legend}} that allow to manipulate the median AUC and legend on the plot
##' @param legend_args character 1L, arguments to \code{\link[base]{text}} and \code{\link[base]{legend}} that allow to manipulate the median AUC and legend on the plot
##' @param ROCR_data logical, return ROC performance data? If FALSE the function only plots ROC curve
##' @return benchmarkMotifsROC: list of ROCR objects (prediction, performance, auc.performance) used to plot ROC curve
##' @import ROCR
##' @export benchmarkMotifsROC
##' @seealso \code{\link{benchmarkMotifs}}
benchmarkMotifsROC = function(res, data_type = c("both","query", "all")[1], col_query = c("grey"), col_all = c("black"), ..., text_args = "col = \"black\"", ROCR_data = F, legend_args = "x = 0.8 | y = 0.25", x_query = 0.25, x_all = 0.75){
  if(data_type == "query") pred = ROCR::prediction(res$predictions_query, res$labels_query)
  if(data_type == "all") pred = ROCR::prediction(res$predictions_all, res$labels_all)
  if(data_type == "both") {
    pred = ROCR::prediction(res$predictions_query, res$labels_query)
    pred_all = ROCR::prediction(res$predictions_all, res$labels_all)
  }

  perf = ROCR::performance(pred, "tpr", "fpr")

  auc.perf = ROCR::performance(pred, measure = "auc")

  if(is.null(text_args)) NULL else {
    text_args = unlist(strsplit(text_args, "\\|"))
  }

  if(data_type == "query"){
    ROCR::plot(perf, col = col_query, ...)
    abline(0,1)
  }
  if(data_type == "all"){
    ROCR::plot(perf, col = col_all, ...)
    abline(0,1)
  }


  if(data_type == "both") {
    ROCR::plot(perf, col = col_query, ...)
    abline(0,1)
    perf_all = ROCR::performance(pred_all, "tpr", "fpr")
    auc.perf_all = ROCR::performance(pred_all, measure = "auc")
    ROCR::plot(perf_all, add = T, col = col_all, ...)
    eval(parse(text = paste0("text(x = x_all, y = median(as.numeric(auc.perf_all@y.values)), labels = paste0(\"Median AUC: \",signif(median(as.numeric(auc.perf_all@y.values)), 3)),",paste0(text_args, collapse = ","),")")))

  }

  eval(parse(text = paste0("text(x = x_query, y = median(as.numeric(auc.perf@y.values)), labels = paste0(\"Median AUC: \",signif(median(as.numeric(auc.perf@y.values)), 3)),",paste0(text_args, collapse = ","),")")))

  if(is.null(legend_args)) NULL else {
    legend_args = unlist(strsplit(legend_args, "\\|"))
  }

  if(data_type == "both"){
    legend_names = c("discovered / known & discoverable",
                     paste0("query instances(proteins): \n",
                            res$N_query_known_instances_found,"(",res$N_query_prot_with_known_instances_found,")",
                            " / ",res$N_query_known_instances,"(",res$N_query_prot_with_known_instances,")"),
                     paste0("non-query instances(proteins): \n",
                            res$N_all_known_instances_found,"(",res$N_all_prot_with_known_instances_found,")",
                            " / ",res$N_all_known_instances,"(",res$N_all_prot_with_known_instances,")"))
    legend_cols = c("transparent",col_query, col_all)
  }
  if(data_type == "all"){
    legend_names = c("discovered / known & discoverable",
                     paste0("non-query instances(proteins): \n",
                            res$N_all_known_instances_found,"(",res$N_all_prot_with_known_instances_found,")",
                            " / ",res$N_all_known_instances,"(",res$N_all_prot_with_known_instances,")"))
    legend_cols = c("transparent", col_all)
  }
  if(data_type == "query"){
  legend_names = c("discovered / known & discoverable",
                   paste0("query instances(proteins): \n",
                          res$N_query_known_instances_found,"(",res$N_query_prot_with_known_instances_found,")",
                          " / ",res$N_query_known_instances,"(",res$N_query_prot_with_known_instances,")"))
  legend_cols = c("transparent",col_query)
  }

  legend_text = paste0("legend(legend = legend_names, col = legend_cols, lty = 1, lwd = 4, merge = TRUE,",
                       paste0(legend_args, collapse = ","), ")")
  eval(parse(text = legend_text))

  if(ROCR_data) return(list(ROC_prediction = pred, ROC_performance = perf, ROC_auc.performance = auc.perf))
}

##' Read QSLIMFinder (SLIMFinder) occurence output to Genomic Ranges object
##' @rdname benchmarkMotifsROC
##' @name mBenchmarkMotifsROC
##' @author Vitalii Kleshchevnikov
##' @param res_list list of objects of class \code{benchmarkMotifsResult}, the output of \code{\link{mBenchmarkMotifs}}
##' @return mBenchmarkMotifsROC: list of lists of ROCR objects (prediction, performance, auc.performance) used to plot ROC curve
##' @export mBenchmarkMotifsROC
##' @seealso \code{\link{benchmarkMotifs}}
mBenchmarkMotifsROC = function(res_list, data_type = "query", ...){
  sapply(res_list, function(res){
    benchmarkMotifsROC(res = res, data_type = data_type, main = res$description, ...)
  })
}
