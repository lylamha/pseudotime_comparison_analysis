library(monocle)
library(data.table)
library(ggplot2)

setwd("/nfs/research2/teichmann/lamha/bifurcation_malaria/monocle")

run_monocle <- function(tpm_file, sample_file, sample_column, out_dir, gene_file = NULL) {
  
  
  # ====================================================
  # get data
  # ====================================================
  
  # expression data
  data <- read.csv(tpm_file, check.names = F, row.names = 1)
  
  # sample annotation
  sampleInfo <- read.csv(sample_file, check.names = F, row.names = 1)
  colnames(sampleInfo)[colnames(sampleInfo) == sample_column] <- "group"
  sampleInfo$group <- as.factor(sampleInfo$group)
  
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
  
  
  # gene annotation
  if (!is.null(gene_file)) {
    geneInfo   <- fread(gene_file)
    geneInfo   <- subset(geneInfo, select=c("Ensembl Gene ID",
                                            "Associated Gene Name",
                                            "Gene type"))
    geneInfo   <- unique(geneInfo)
    
    
    tmp <- data.frame(genes = rownames(data), seq = 1:nrow(data))
    geneInfo <- merge(geneInfo, tmp ,
                      by.x="Ensembl Gene ID", by.y="genes", all.y=T)
    geneInfo <- geneInfo[order(geneInfo$seq),]
    geneInfo$seq <- NULL
    geneInfo <- as.data.frame(geneInfo)
    rownames(geneInfo) <- geneInfo$`Ensembl Gene ID`
    
  } else {
    geneInfo <- data.frame("Ensembl Gene ID" = rownames(data),
                           check.names = F)
    rownames(geneInfo) <- rownames(data)
  }
  
  
  fd <- new("AnnotatedDataFrame", data = geneInfo)
  
  
  # ====================================================
  # create cell data set and transform
  # ====================================================
  
  
  pd <- new("AnnotatedDataFrame", data = sampleInfo)
  
  cell_data <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 1, 
                              expressionFamily = tobit())
 
  
  cell_data <- estimateSizeFactors(cell_data)
  # cell_data <- estimateDispersions(cell_data)
  
  
  # detect genes
  cell_data <- detectGenes(cell_data, min_expr = 0.1)
  
  
  # filter by number of expressed genes
  expressed_genes <- row.names(subset(fData(cell_data),
                                      num_cells_expressed >= 50))
  
  
  # log transform expression data
  L <- log(exprs(cell_data[expressed_genes,]))
  
  
  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily"
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  
  
  # Plot the distribution of the standardized gene expression values.
  qplot(value, geom="density", data=melted_dens_df) +
    stat_function(fun = dnorm, size=0.5, color="red") +
    xlab("Standardized log(FPKM)") +
    ylab("Density")
  ggsave(paste0(out_dir, "/gene_expression_log_distribution.pdf"))
  
  
  # ====================================================
  # differential expression 
  # ====================================================
  
  diff_test_res <- differentialGeneTest(cell_data[expressed_genes,],
                                        fullModelFormulaStr="~group")
  
  ### significant genes
  
  sig_genes <- subset(diff_test_res, qval < 0.1)
  head(sig_genes[order(sig_genes$qval),], n=20)
  
  # ====================================================
  # ordering cells by progress
  # ====================================================
  
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  ordering_genes <- intersect(ordering_genes, expressed_genes)
  
  # filter cells
  cell_data <- setOrderingFilter(cell_data, ordering_genes)
  
  # reduce dimension
  cell_data <- reduceDimension(cell_data, max_components=2)
  
  # order cells
  cell_data <- orderCells(cell_data, reverse=FALSE)
  
  
  # ====================================================
  # plot trajectory
  # ====================================================
  
  plot_cell_trajectory(cell_data, color_by = "group") +
  ggsave(paste0(out_dir, "/cell_trajectory_sampleInfo.pdf"))
  
  
  plot_cell_trajectory(cell_data, color_by = "State")
  ggsave(paste0(out_dir, "/cell_trajectory_stateInfo.pdf"))
  
  
}


# ======================================
# run monocle
# ======================================

run_monocle("../data/malaria/all_tpm_QCfiltered.csv",
            "../data/malaria/tcells_rebuttal.csv", 
            "day",
            "malaria")

run_monocle("../data/lung/lung_tpm.csv",
            "../data/lung/lung_sample_info.csv",
            "condition",
            "lung")

run_monocle("../data/frog/frog_tpm.csv",
            "../data/frog/frog_sample_info.csv",
            "hpf",
            "frog")

run_monocle("../data/pgc/pgc_tpm.csv",
            "../data/pgc/pgc_sample_info.csv",
            "week",
            "pgc")

