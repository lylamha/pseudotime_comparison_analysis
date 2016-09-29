
library(Mpath)

### preprocessing tpm files (remove ERCCs and rescale)
### generate tab-delimited data and sample annotation files
### sample info: 2 columns (cell, GroupID)
### returns new tpm and sample filenames
data_init <- function(tpm_file, sample_file, sample_column) {
  
  new_tpm_file    <- gsub(".csv", "_mpath.txt", tpm_file)
  new_sample_file <- gsub(".csv", "_mpath.txt", sample_file)
  
  if (! (file.exists(new_tpm_file) && file.exists(new_sample_file)) ) {
    
    # expression data
    data <- read.csv(tpm_file, check.names = F, row.names = 1)
    
    # sample annotation
    sampleInfo <- read.csv(sample_file, check.names = F, row.names = 1)
    colnames(sampleInfo)[colnames(sampleInfo) == sample_column] <- "GroupID"
    sampleInfo$cell <- rownames(sampleInfo)
    sampleInfo <- subset(sampleInfo, select = c(cell, GroupID))
    
    # remove NA's
    rm.idx <- which(colSums(apply(data, 2, is.na)) > 0)
    
    if (length(rm.idx)) {
      data <- data[,-rm.idx]
      sampleInfo <- sampleInfo[-rm.idx,]
    }
    
    ### exclude ERCC and rescale
    data    <- data[-grep("ERCC", rownames(data)),]
    rescale <- colSums(data)
    data <- t(data) / rescale * 1e6
    data <- t(data)
    
    # write out data for mpath
    print("generate new tpm and sample file")
    write.table(data, new_tpm_file, sep="\t", row.names = T, col.names = T, quote = F)
    
    write.table(sampleInfo, new_sample_file, sep="\t", row.names = F, col.names = T, quote = F)
    
  }
  
  return(list(tpm_file = new_tpm_file, sample_file = new_sample_file))
  
}


run_mpath <- function(tpm_file, sample_file, sample_column, out_dir,
                      optimal_number_cluster = NULL) {

  meta_dir <- paste0(out_dir, "/mpath_metadata")
  if (! file.exists( meta_dir )) { dir.create(meta_dir) }
  
  # =============================================
  # preprocess data for Mpath
  # =============================================
  
  data_files <- data_init(tpm_file, sample_file, sample_column)
  
  # =============================================
  # quality check
  # =============================================
  
  # QC_gene removes genes that have TPM values < 1
  # in more than 95 percent of cells in each group
  
  rpkmFile   <- data_files$tpm_file
  file_bname <- gsub(".txt", "", basename(data_files$tpm_file))
  rpkmQCFile <- paste0(meta_dir, "/", file_bname, "_geneQC.txt")
  
  if (out_dir == "frog") {
    sFile <- NULL
  } else {
    sFile <- data_files$sample_file
  }
  QC_gene(rpkmFile   = rpkmFile,
          rpkmQCFile = rpkmQCFile,
          sampleFile = sFile,
          threshold  = 0.05, method = "any")
  
  sampleFile <- data_files$sample_file
  
  # =============================================
  # find optimal cluster
  # =============================================
  
  if (is.null(optimal_number_cluster)) {
    plot_fname <- sub(".txt", "", rpkmQCFile)
    plot_fname <- paste0(plot_fname, "_ncluster_vs_nlm.pdf")
    # generates plots where rpkmQCFile is located
    find_optimal_cluster_number(rpkmFile = rpkmQCFile,
                                sampleFile = sampleFile,
                                min_cluster_num = 7,
                                max_cluster_num = 20,
                                diversity_cut = 0.6,
                                size_cut = 0.05)
    return(paste0("please inspect ", 
                  plot_fname,
                  " to get optimal number of cluster",
                  " and to continue"))
    
  }
  
  # =============================================
  # Landmark designation
  # =============================================
  
  baseName   <- sub(".txt", "", rpkmQCFile)
  landmark_cluster <- landmark_designation(rpkmFile = rpkmQCFile,
                                           baseName = baseName,
                                           sampleFile = sampleFile,
                                           method = "diversity_size",
                                           numcluster = optimal_number_cluster,
                                           diversity_cut=0.6,
                                           size_cut=0.05)
  
  ### Plot hierachical clustering 
  SC_hc_colorCode(dataFile   = rpkmQCFile,
                  cuttree_k  = optimal_number_cluster,
                  sampleFile = sampleFile,
                  width = 22, height = 10, iflog2 = TRUE,
                  colorPalette = NULL)
  
  
  # =============================================
  # Construct weighted neighborhood network
  # =============================================

  neighbor_network <- build_network(exprs = rpkmQCFile,
                                    landmark_cluster = landmark_cluster,
                                    baseName = baseName)
  
  
  ### TrimNet: trim edges of lower weights
  trimmed_net <- trim_net(neighbor_network, textSize=30,
                          baseName = baseName,
                          method = "mst")
  
  
  
  return(1)

}

run_mpath(tpm_file      = "../data/lung/lung_tpm.csv",
          sample_file   = "../data/lung/lung_sample_info.csv",
          sample_column = "condition",
          out_dir       = "lung",
          optimal_number_cluster = 19)


run_mpath(tpm_file      = "../data/frog/frog_tpm.csv",
          sample_file   = "../data/frog/frog_sample_info.csv",
          sample_column = "hpf",
          out_dir       = "frog",
          optimal_number_cluster = NULL)


run_mpath(tpm_file      = "../data/pgc/pgc_tpm.csv",
          sample_file   = "../data/pgc/pgc_sample_info.csv",
          sample_column = "week",
          out_dir       = "pgc",
          optimal_number_cluster = 24)


run_mpath(tpm_file      = "../data/malaria/all_tpm_QCfiltered.csv",
          sample_file   = "../data/malaria/tcells_rebuttal.csv",
          sample_column = "day",
          out_dir       = "malaria",
          optimal_number_cluster = 10)

