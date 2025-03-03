library(scSHC)
type <- "RNA"
method <- "scSHC-RNA"

if (dim(temp_object@assays$RNA@counts)[1] >= 2500) {  
  temp_num_features <- 2500  
} else {  
  temp_num_features <- dim(temp_object@assays$RNA@counts)[1]  
} 

start_time <- Sys.time()
mem <- peakRAM({
  scSHC_cluster_result <- scSHC(temp_object@assays$RNA@counts, num_features=temp_num_features)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
pred_label <- scSHC_cluster_result[[1]]
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

detach("package:scSHC", unload=TRUE)
