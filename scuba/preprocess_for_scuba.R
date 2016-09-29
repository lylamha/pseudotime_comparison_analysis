
preprocess_for_SCUBA <- function(tpm_file, sample_file, sample_column) {
  
  data <- read.csv(tpm_file, row.names = 1, check.names=F, stringsAsFactors = F)
  
  sampleInfo <- read.csv(sample_file, check.names=F, stringsAsFactors = F)
  
  ### remove NA's
  rm.idx <- which(colSums(apply(data, 2, is.na)) > 0)
  if (length(rm.idx)) {
    data       <- data[,-rm.idx]
    sampleInfo <- sampleInfo[-rm.idx,]
  }
  
  ### exclude ERCC and rescale
  data    <- data[-grep("ERCC", rownames(data)),]
  rescale <- colSums(data)
  data <- t(data) / rescale * 1e6
  data <- t(data)
  
  data <- cbind("Cell ID" = rownames(data),
                data)
  
  if (length(grep("lung", tpm_file))) {
    sampleInfo[,sample_column] <- factor(sampleInfo[,sample_column], level=c("E14.5",
                                                      "E16.5",
                                                      "E18.5",
                                                      "AT2"))
    
  }
  
  newrow <- c("stage", sampleInfo[,sample_column])
  
  data <- rbind(newrow, data)
  
  file <- gsub(".csv", "_scuba.txt", tpm_file)
  write.table(data, file, sep="\t", row.names = F, col.names = T, quote = F)
}



preprocess_for_SCUBA("../data/malaria/all_tpm_QCfiltered.csv",
                     "../data/malaria/tcells_rebuttal.csv",
                     "day")

preprocess_for_SCUBA("../data/lung/lung_tpm.csv",
                     "../data/lung/lung_sample_info.csv",
                     "condition")

preprocess_for_SCUBA("../data/frog/frog_tpm.csv",
                     "../data/frog/frog_sample_info.csv",
                     "hpf")

preprocess_for_SCUBA("../data/pgc/pgc_tpm.csv",
                     "../data/pgc/pgc_sample_info.csv",
                     "week")
