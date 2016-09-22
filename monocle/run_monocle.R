library(monocle)
library(data.table)
library(ggplot2)

data <- read.delim("monocle/only_genes_CD4tcells_tpm.txt", check.names = F, row.names = 1)
sampleInfo <- read.delim("monocle/CD4tcells_sample_day_info.txt", check.names = F, row.names = 1)
geneInfo   <- fread("/nfs/research2/teichmann/reference/mus-musculus/ensembl85-mouse-tx-annotation.csv")
geneInfo   <- subset(geneInfo, select=c("Ensembl Gene ID",
                                        "Associated Gene Name",
                                        "Gene type"))
geneInfo   <- unique(geneInfo)

# gene annotation

tmp <- data.frame(genes = rownames(data), seq = 1:nrow(data))
geneInfo <- merge(geneInfo, tmp , by.x="Ensembl Gene ID", by.y="genes", all.y=T)
geneInfo <- geneInfo[order(geneInfo$seq),]
geneInfo$seq <- NULL
geneInfo <- as.data.frame(geneInfo)
rownames(geneInfo) <- geneInfo$`Ensembl Gene ID`


# create cell data set

pd <- new("AnnotatedDataFrame", data = sampleInfo)
fd <- new("AnnotatedDataFrame", data = geneInfo)
cell_data <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit=1, 
                            expressionFamily = tobit())

cell_data <- estimateSizeFactors(cell_data)
# cell_data <- estimateDispersions(cell_data)


# detect genes

cell_data <- detectGenes(cell_data, min_expr = 0.1)

# filter by number of expressed genes

expressed_genes <- row.names(subset(fData(cell_data), num_cells_expressed >= 50))

L <- log(exprs(cell_data[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color="red") +
  xlab("Standardized log(FPKM)") +
  ylab("Density")
ggsave(paste0(out_dir, "/gene_expression_log_distribution.pdf"))

