library(cidr)
type <- "ADT"
method <- "CIDR-ADT"

start_time <- Sys.time()
mem <- peakRAM({
  scCIDR <- scDataConstructor(as.matrix(temp_object@assays$ADT@counts))
  scCIDR <- determineDropoutCandidates(scCIDR)
  scCIDR <- wThreshold(scCIDR)
  scCIDR <- scDissim(scCIDR)
  scCIDR <- scPCA(scCIDR)
  scCIDR <- nPC(scCIDR)
  scCIDR <- scCluster(scCIDR,nCluster=n_cluster)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
pred_label <- scCIDR@clusters
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

detach("package:cidr", unload=TRUE)