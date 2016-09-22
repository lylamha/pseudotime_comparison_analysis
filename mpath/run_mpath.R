

library(Mpath)

# =============================================
# quality check
# =============================================

# QC_gene removes genes that have TPM values < 1
# in more than 95 percent of cells in each group

rpkmFile   <- "mpath/only_genes_CD4tcells_tpm.txt"
rpkmQCFile <- "mpath/only_genes_CD4tcells_tpm_geneQC.txt"
sampleFile <- NULL
QC_gene(rpkmFile   = rpkmFile,
        rpkmQCFile = rpkmQCFile,
        sampleFile = sampleFile, threshold = 0.05, method = "any")