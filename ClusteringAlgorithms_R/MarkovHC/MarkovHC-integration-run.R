library(MarkovHC)
library(stringr)
type <- "integration_feature"
method <- "MarkovHC-integration_feature"

temp_dim <- dim(temp_object@reductions$pca_adt)[2]
temp_dim <- ifelse(temp_dim < 9, temp_dim, 9)
MarkovHC_temp_object <- FindNeighbors(object = temp_object,
                                      k.param = 70,
                                      compute.SNN = TRUE,
                                      prune.SNN = 0,
                                      reduction = "pca_adt",
                                      dims = 1:temp_dim,
                                      force.recalc = TRUE)

start_time <- Sys.time()
mem <- peakRAM({
  MarkovHC_object <- MarkovHC(MarkovHC_input = MarkovHC_temp_object,
                              dobasecluster = TRUE,
                              SNNslot = 'ADT_snn',
                              KNNslot = 'ADT_nn',
                              cutpoint = 0.001,
                              verbose = FALSE)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")


internal_measures <- IMI_selection(MarkovObject=MarkovHC_object,
                                   prune=TRUE,
                                   weed=10)
MarkovHC_culster_result <-  fetchLabels(MarkovObject=MarkovHC_object,
                                        MarkovLevels=1:length(MarkovHC_object$hierarchicalStructure),
                                        prune = TRUE, weed = 10)


true_label <- factor(get(temp_dataset)@meta.data$celltype)

#####################################################################################################
##Select levels that match the number of true cell types.
column_names <- names(MarkovHC_culster_result)
pred_label <- NULL

for (i in 1:length(column_names)) {
  if(length(unique(MarkovHC_culster_result[, i]))==n_cluster)
    pred_label <- MarkovHC_culster_result[, i]
}

if (is.null(pred_label)) {
  print("Unable to find results that meet the criteria!!!")
  temp_metric <- data.frame(ari=-1, nmi=-1, purity=-1, acc=-1)
} else {
  pred_label_for_metric_compute <- pred_label
  values_with_plus <- grep("\\+", pred_label_for_metric_compute, value = TRUE)
  values <- as.numeric(pred_label_for_metric_compute[!grepl("\\+", pred_label_for_metric_compute)])
  max_value <- max(values, na.rm = TRUE)
  new_values <- seq(max_value + 1, length.out = length(values_with_plus))
  for (i in seq_along(values_with_plus)) {
    pred_label_for_metric_compute[pred_label_for_metric_compute == values_with_plus[i]] <- new_values[i]
  }
  temp_metric <- metric_compute(true_label, pred_label_for_metric_compute)
}
#####################################################################################################

file_name <- paste(method, "_pred.csv", sep = "")
pred_csv <- file.path(result_dir, file_name)
result_data <- data.frame(pred_label)
write.csv(result_data, file = pred_csv, row.names = FALSE)

result_row <- data.frame(t(c(dataset_, type, method, temp_metric$ari, temp_metric$nmi, temp_metric$purity, temp_metric$acc, peak_mem, running_time_temp)))

write.table(
  result_row,
  file = result_csv,
  append = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  sep = ","
)

detach("package:MarkovHC", unload=TRUE)
detach("package:stringr", unload=TRUE)