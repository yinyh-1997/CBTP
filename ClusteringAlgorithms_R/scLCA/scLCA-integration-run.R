library(scLCA)
type <- "integration_feature"
method <- "scLCA-integration_feature"

if (dim(temp_object@assays$ADT@counts)[1] == 19) {
  temp_sclca_input <- rbind(as.matrix(temp_object@assays$ADT@counts),as.matrix(temp_object@assays$ADT@counts))
} else {
  temp_sclca_input <- as.matrix(temp_object@assays$ADT@counts)
}

start_time <- Sys.time()
mem <- peakRAM({
  scLCA_cluster_result <- myscLCA(temp_sclca_input, clust.max = n_cluster+5)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
true_label <- factor(temp_object@meta.data$celltype)
pred_label <- scLCA_cluster_result[[1]]
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

detach("package:scLCA", unload=TRUE)
