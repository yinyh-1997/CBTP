library(SC3)
library(SingleCellExperiment)
library(scater)
type <- "integration_feature"
method <- "SC3-integration_feature"

start_time <- Sys.time()
mem <- peakRAM({
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(temp_object@assays$ADT@counts),
      logcounts = as.matrix(temp_object@assays$ADT@counts)
    )
  )
  
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  sce <- runPCA(sce)
  
  sce <- sc3(sce,ks = n_cluster:n_cluster,svm_max=dim(temp_object@assays$ADT@counts)[2],rand_seed=benchmark_seed,gene_filter = FALSE)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

all_column_names <- names(sce@colData@listData)
true_label <- factor(get(temp_dataset)@meta.data$celltype)
for (col_name in all_column_names) {
  pred_label <- sce@colData@listData[[col_name]]
}
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
