drugs <- fread("/data1/meaneylab/eamon/Placenta_WGCNA/Drug_analysis/Drug_gene_DGdb_April_29th_2024")
DGI <- split(drugs$gene_name,drugs$drug_name)


load("/data1/meaneylab/eamon/Placenta_WGCNA/rWGCNA/Reup/Two_batches/Two_batches_filt_norm_corr.Rdata")

ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

lookup <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "external_gene_name"),
                values=rownames(datExpr.filt.norm.corr),mart= ensembl_human) 

# Ensure rownames are characters for accurate matching
rownames(datExpr.filt.norm.corr) <- as.character(rownames(datExpr.filt.norm.corr))

# Use `match` to find the indices of the old rownames in lookup_df
indices <- match(rownames(datExpr.filt.norm.corr), lookup$ensembl_gene_id)

# Replace rownames with the new names from lookup_df
rownames(datExpr.filt.norm.corr) <- lookup$external_gene_name[indices]

library(GSVA)
library(BiocParallel)
# Register BiocParallel to use multiple cores
bpparam <- MulticoreParam(workers = 5)  # Adjust the number of workers based on your system’s capabilities

# Compute ssGSEA scores
ssgsea_scores <- gsva(datExpr.filt.norm.corr, DGI,
                      method = "ssgsea",
                      min.sz = 10,
                      max.sz = 500,
                      BPPARAM = bpparam)

ssgsea_scores <- as.data.frame(t(ssgsea_scores))
ssgsea_scores <- rownames_to_column(ssgsea_scores, var = "Sample")

# PCs
load("/data1/meaneylab/eamon/Placenta_WGCNA/ML/NDD_PC_generation.RData")
PCs <- p$rotated
variance_exp <- cumsum(p$variance)

# 82PCs explain 98.8% of variance
PCs <- PCs[,1:82]
PCs <- rownames_to_column(PCs, var = "Sample")


# Merge the dataframes by their row names
combined_data <- merge(PCs, ssgsea_scores, by = "Sample", all = TRUE)
combined_data$Sample <- NULL

# Correlation matrix
# Assuming PCs columns start with "PC" and ssgsea_scores columns with "GSE"
pc_columns <- grep("^PC", names(combined_data), value = TRUE)
ssGSEA_names <- setdiff(names(ssgsea_scores), "Sample")

# List to store dataframes
list_of_dfs <- list()

# Total number of PCs
total_pcs <- length(pc_columns)

for (index in seq_along(pc_columns)) {
  pc <- pc_columns[index]
  
  # Print progress update
  cat(sprintf("Processing %s (%d of %d)\n", pc, index, total_pcs))
  
  # Prepare a dataframe to hold results for this PC
  results_df <- data.frame(
    Term = ssGSEA_names,
    Correlation = NA_real_,
    P_value = NA_real_,
    P_value_Bonferroni = NA_real_,
    P_value_FDR = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Store raw P-values temporarily
  raw_p_values <- numeric(length(ssGSEA_names))
  
  # Calculate correlation and P-value for each GSE score
  for (i in seq_along(ssGSEA_names)) {
    gse <- ssGSEA_names[i]
    test_result <- cor.test(combined_data[[pc]], combined_data[[gse]], 
                            method = "pearson", use = "complete.obs")
    
    # Store results in the dataframe
    results_df$Correlation[i] <- test_result$estimate
    raw_p_values[i] <- test_result$p.value
  }
  
  # Adjust P-values
  results_df$P_value <- raw_p_values
  results_df$P_value_Bonferroni <- p.adjust(raw_p_values, method = "bonferroni")
  results_df$P_value_FDR <- p.adjust(raw_p_values, method = "fdr")
  
  # Add the results dataframe to the list with the PC name as the list key
  list_of_dfs[[pc]] <- results_df
}

Check_me <- list_of_dfs$PC20
