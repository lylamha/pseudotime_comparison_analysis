
library(dpt)
library(ggplot2)

setwd("/nfs/research2/teichmann/lamha/bifurcation_malaria/dpt")

run_dpt <- function(tpm_file, sample_file, sample_column, root_cell, out_dir) {
  
  # =============================================================
  # TPM data
  # =============================================================
  
  data <- read.csv(tpm_file, row.names=1, check.names=F)
  data <- t(data) ## transform to cells x genes
  
  ### exclude ERCC and rescale
  data    <- data[,-grep("ERCC", colnames(data))]
  rescale <- rowSums(data)
  data <- data / rescale * 1e6
  
  # =============================================================
  # sample information
  # =============================================================
  
  sampleInfo <- read.csv(sample_file, check.names=F, row.names = 1)
  
  ### set colors ?
  # mycolors <- c("#81CDC3", "#EB867F",  "#D0AB80", "#BDD68F", "#F9EF76")
  # names(mycolors) <- c("day_0", "day_2", "day_3", "day_4", "day_7")
  # 
  # colScale <- scale_colour_manual(name = "days", values = mycolors)
  
  # =============================================================
  # finding sigma for gaussian kernel
  # =============================================================
  
  sigmas <- destiny::find.sigmas(data, verbose = F)
  
  # =============================================================
  # execute DPT
  # =============================================================
  
  ts <- Transitions(data, sigma = sigmas@optimal.sigma)
  
  root_idx <- which(rownames(data) == root_cell)
  
  
  # Diffusion pseudotime
  pt <- dpt(ts, branching = TRUE, root = root_idx)
  
  
  ev <- eigen(as.matrix(ts@transitions), TRUE)$vectors
  dm <- as.data.frame(ev[, -1])
  colnames(dm) <- paste0('DC', seq_len(ncol(dm)))
  
  
  # =============================================================
  # output
  # =============================================================
  
  qplot(DC1, DC2, data = dm, colour = pt$Branch) +
    theme_bw()
  ggsave(paste0(out_dir, "/branching_tree_branchInfo.pdf"))
  
  
  qplot(DC1, DC2, data = dm, colour = as.factor(sampleInfo[,sample_column])) +
    scale_fill_continuous(guide = guide_legend(title = sample_column)) +
    theme_bw()
  ggsave(paste0(out_dir, "/branching_tree_sampleInfo.pdf"))
  
  return(1)
}

# ======================================
# run dpt
# ======================================

run_dpt("../data/lung/lung_tpm.csv",
        "../data/lung/lung_sample_info.csv",
        "condition",
        "SRP033209_E14.5_rep_1_cell_24",
        "lung")

run_dpt("../data/frog/frog_tpm.csv",
        "../data/frog/frog_sample_info.csv",
        "hpf",
        "SRR1795631",
        "frog")

run_dpt("../data/pgc/pgc_tpm.csv",
        "../data/pgc/pgc_sample_info.csv",
        "week",
        "SRR2013600",
        "pgc")




