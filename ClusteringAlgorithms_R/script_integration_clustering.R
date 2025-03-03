######################################################################################################################
######################################Definition of Calculation Functions#############################################
######################################################################################################################
# Load necessary libraries for linear matching, ARI, and NMI calculations, memory peak tracking, and Seurat.
library(clue)
library(aricode)
library(peakRAM)
library(Seurat)

benchmark_cluster_acc <- function(y_true, y_pred) {
  # Ensure y_true and y_pred are integer vectors
  y_true <- as.integer(y_true)
  y_pred <- as.integer(y_pred)
  
  # Check that y_true and y_pred have the same length
  if (length(y_true) != length(y_pred)) {
    stop("y_true and y_pred must have same length!")
  }
  
  # Calculate the maximum value in y_true and y_pred
  D <- max(max(y_true), max(y_pred)) + 1

  # Initialize a matrix with zeros
  w <- matrix(0, nrow = D, ncol = D)
  
  # Count occurrences
  for (i in 1:length(y_true)) {
    w[y_true[i] + 1, y_pred[i] + 1] <- w[y_true[i] + 1, y_pred[i] + 1] + 1
  }
  
  # Perform a linear assignment problem
  ind <- clue::solve_LSAP(t(w), maximum = TRUE)
  
  # Return the accuracy
  return(sum(w[cbind(ind, 1:length(ind))]) / length(y_pred))
}
#########################################
#########################################
benchmark_purity <- function(y_true, y_pred){
  
  # Convert to vectors
  y_true <- as.vector(y_true)
  y_pred <- as.vector(y_pred)
  
  # Retrieve unique predicted clusters
  clusters <- unique(y_pred)
  n_clusters <- length(clusters)
  
  # Initialize storage for each cluster's labels
  cluster_labels <- vector(length=n_clusters) 
  
  # Calculate the most common true label for each cluster
  for(i in seq_along(clusters)){
    cluster_mask <- y_pred == clusters[i]
    if(sum(cluster_mask) > 0){
      cluster_labels[i] <- names(which.max(table(y_true[cluster_mask])))  
    }
  }
  
  # Count correct classifications
  correct <- 0
  for(i in seq_along(y_true)){
    if(!is.na(match(y_pred[i], clusters))){
      correct <- correct + as.integer(y_true[i] == cluster_labels[match(y_pred[i], clusters)])
    } 
  }
  
  # Calculate purity
  purity <- correct / length(y_true)
  
  return(purity)
}

metric_compute <- function(y_true, y_pred) {
  temp_ari <- aricode::ARI(y_true, y_pred)
  temp_nmi <- aricode::NMI(y_true, y_pred,variant="sum")
  # The use of variant="sum" indicates that the arithmetic mean normalization method is employed. 
  # This is consistent with the average_method="arithmetic" in the Python function normalized_mutual_info_score. 
  # Therefore, in this case, the NMI computation results from R and Python should be equivalent.
  temp_purity <- benchmark_purity(y_true, y_pred)
  temp_acc <- benchmark_cluster_acc(y_true, y_pred)
  metrics <- data.frame(ari=temp_ari, nmi=temp_nmi, purity=temp_purity, acc=temp_acc)
  return(metrics)
}
######################################################################################################################
##############################################Definition of Paths#####################################################
######################################################################################################################
dataset_path <- "/home/rstudio/public_data/AllDataset/Dataset_integration/run" 
current_time <- Sys.time()
current_time_str <- gsub(" ","_",gsub(" ","_",format(current_time, "%Y-%m-%d %H-%M-%S")))
dataset_list <- list.dirs(dataset_path, full.names = FALSE, recursive = FALSE)
result_path <- "/home/rstudio/public_data/BenchmarkResult_integration"
benchmark_seed <- 2024
######################################################################################################################
################################Loop Through Datasets to Execute Various Algorithms###################################
######################################################################################################################

for (dataset_ in dataset_list) {

  print(dataset_)
  # Concatenate the full dataset path
  dataset <- file.path(dataset_path, dataset_)
  temp_dataset <- paste0(dataset_, "_object")
  assign(x = temp_dataset, value = readRDS(paste0(dataset, "/CITE.rds")))
  # Concatenate the result saving path
  result_dir <- file.path(result_path, current_time_str, dataset_)
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }

  result_csv <- file.path(result_dir, "Result.csv")
  # Write the header if the CSV does not exist
  if (!file.exists(result_csv)) {
    write.table(
      data.frame(Dataset = character(0), Type = character(0), Method = character(0), ARI = character(0), NMI = character(0), Purity = character(0), Cluster_ACC = character(0), Memory_peak = character(0), Running_time = character(0)),
      file = result_csv,
      row.names = FALSE,
      col.names = c("Dataset", "Type", "Method", "ARI", "NMI", "Purity", "ACC", "Memory_peak", "Running_time"),
      sep = ","
    )
  }

  temp_object <- get(temp_dataset)
  var_genes_matrix <- GetAssayData(temp_object, layer = "data")
  n_cluster <- length(unique(temp_object@meta.data$celltype))
  
  
  # ###############################################################################################
  # Load R script files

  ###Celda
  local({tryCatch({
    source("/home/rstudio/public_data/Celda/Celda-integration-run.R", local = TRUE)  }, error = function(e) {} )})
  
  ###CIDR
  local({tryCatch({
    source("/home/rstudio/public_data/CIDR/CIDR-integration-run.R", local = TRUE)  }, error = function(e) {} )})
  
  ###scSHC
  local({tryCatch({
    source("/home/rstudio/public_data/scSHC/scSHC-integration-run.R", local = TRUE)  }, error = function(e) {} )})
  
  ###scLCA
  local({tryCatch({
    source("/home/rstudio/public_data/scLCA/scLCA-integration-run.R", local = TRUE) }, error = function(e) {} )})
  
  ###Monocle3
  local({tryCatch({
    source("/home/rstudio/public_data/Monocle3/Monocle3-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###FlowSOM
  local({tryCatch({
    source("/home/rstudio/public_data/FlowSOM/FlowSOM-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###TSCAN
  local({tryCatch({
    source("/home/rstudio/public_data/TSCAN/TSCAN-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###DEPECHER
  local({tryCatch({
    source("/home/rstudio/public_data/DEPECHER/DEPECHER-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###SHARP
  local({tryCatch({
    source("/home/rstudio/public_data/SHARP/SHARP-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###DRSC
  local({tryCatch({
    source("/home/rstudio/public_data/DRSC/DRSC-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###MarkovHC
  local({tryCatch({
    source("/home/rstudio/public_data/MarkovHC/MarkovHC-integration-run.R", local = TRUE)  }, error = function(e) {} )})
  
  ###CDC
  local({tryCatch({
    source("/home/rstudio/public_data/CDC/CDC-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###SIMLR
  local({tryCatch({
    source("/home/rstudio/public_data/SIMLR/SIMLR_Large_Scale-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###Spectrum
  local({tryCatch({
    source("/home/rstudio/public_data/Spectrum/Spectrum-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ###SC3
  local({tryCatch({
    source("/home/rstudio/public_data/SC3/SC3-integration-run.R", local = TRUE)  }, error = function(e) {} )})

  ## After each iteration, remove the currently loaded dataset to free memory
  rm(list = temp_dataset)
}
