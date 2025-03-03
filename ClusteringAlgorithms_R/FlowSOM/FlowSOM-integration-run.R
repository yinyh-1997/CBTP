library(FlowSOM)
type <- "integration_feature"
method <- "FlowSOM-integration_feature"

start_time <- Sys.time()
mem <- peakRAM({
  FlowSOM_cluster_result <- FlowSOM(t(as.matrix(temp_object@assays$ADT@data)),nClus=n_cluster, seed=benchmark_seed)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
pred_label <- GetMetaclusters(FlowSOM_cluster_result)
temp_metric <- metric_compute(true_label, pred_label)

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

detach("package:FlowSOM", unload=TRUE)

