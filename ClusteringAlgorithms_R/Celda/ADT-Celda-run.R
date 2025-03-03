library(celda)
type <- "ADT"
method <- "Celda-ADT"


temp_counts_matrix <- temp_object@assays$ADT@counts
temp_non_zero_adts <- rowSums(temp_counts_matrix != 0) > 0
DefaultAssay(temp_object) <- "ADT" 
temp_object <- subset(temp_object, features = rownames(temp_counts_matrix)[temp_non_zero_adts])

start_time <- Sys.time()
mem <- peakRAM({
  Celda_cluster_result <- celda_C(x = temp_object@assays$ADT@counts, K = n_cluster, seed = benchmark_seed)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
pred_label<-celdaClusters(Celda_cluster_result)
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

detach("package:celda", unload=TRUE)
