library('DR.SC')
type <- "integration_feature"
method <- "DR.SC-integration_feature"

temp_row <- dim(temp_object@assays$ADT@data)[1]

if (temp_row > 15) {
  temp_input <- t(temp_object@assays$ADT@data)
  temp_q <- 15
} else if (temp_row == 10) {
  temp_input <- t(temp_object@assays$ADT@data)
  temp_q <- 9
} else if (temp_row == 2) {
  temp_input <- t(rbind(temp_object@assays$ADT@data, temp_object@assays$ADT@data))
  temp_q <- 2
}

start_time <- Sys.time()
mem <- peakRAM({
  DRSC_result <- DR.SC_fit(X=temp_input, K=n_cluster,q=temp_q)
})
peak_mem <- mem$Peak_RAM_Used_MiB
end_time <- Sys.time()
running_time_temp <- as.numeric(end_time - start_time, units = "secs")

true_label <- factor(get(temp_dataset)@meta.data$celltype)
pred_label <- c(DRSC_result$Objdrsc_nonspa[[1]]$cluster)
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

detach("package:DR.SC", unload=TRUE)

