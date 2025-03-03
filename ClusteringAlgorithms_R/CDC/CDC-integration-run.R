local({

type <- "integration_feature"
method <- "CDC-integration_feature"
source('/home/rstudio/public_data/CDC.R', local=TRUE)
k = 30 
ratio = 0.9 

temp_dim <- dim(temp_object@reductions$pca_adt)[2]
temp_dim <- ifelse(temp_dim < 20, temp_dim, 20)
temp_CDC_ADT_object <- RunUMAP(temp_object, reduction = "pca_adt", dims = 1:temp_dim, assay = "ADT")

start_time <- Sys.time()
mem <- peakRAM({
  CDC_cluster_result <- CDC(temp_CDC_ADT_object@reductions$umap@cell.embeddings, k, ratio)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
pred_label <- as.factor(CDC_cluster_result)
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

})
