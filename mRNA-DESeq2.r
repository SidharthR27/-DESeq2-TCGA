#Intial Setup

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
BiocManager::install("DESeq2")
library(DESeq2)

#Function for Identifiying Upregulated and Downregulated genes

get_upregulated <- function(df){
    key <- intersect(rownames(df)[which(df$log2FoldChange>=2)],
              rownames(df)[which(df$pvalue==0.05)])
    
    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    return(results)
}
get_downregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-2)],
            rownames(df)[which(df$pvalue==0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

#Downloading and Data Wrangling

query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
GDCdownload(query, method = "api", files.per.chunk = 100,
            directory = "C:/Users/thero/Documents/content/Data")
mrna_df <- GDCprepare(query, directory = "C:/Users/thero/Documents/content/Data")

# Remove columns we dont need, keep counts

mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_df <- assay(mrna_df)
delim_fn = function(x, n, i){
    do.call(c, lapply(x, function(X)
        paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}
colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 4)
mrna_meta <- as.data.frame(mrna_meta)
mrna_df <- as.data.frame(mrna_df)

# Remove metastatic sample

metastatic_key <- mrna_meta[which(mrna_meta[,2] == "Metastatic"),]
mrna_meta <- mrna_meta[!mrna_meta[,2] == metastatic_key[,2],]
mrna_df <- mrna_df[, -grep(paste0(metastatic_key[,1]), colnames(mrna_df))]
mrna_meta[,2] <- as.character(mrna_meta[,2])
mrna_meta[,2] <- gsub("Primary solid Tumor", "Tumor", mrna_meta[,2])
mrna_meta[,2] <- gsub("Solid Tissue Normal", "Normal", mrna_meta[,2])
mrna_meta[,2] <- as.factor(mrna_meta[,2])
levels(mrna_meta[,2])
colnames(mrna_meta) <- c("cases", "Condition")

# DESeq2 Analysis

mrna_dds <- DESeqDataSetFromMatrix(round(mrna_df), colData = mrna_meta, design = ~ Condition)
mrna_dds$Condition <- relevel(mrna_dds$Condition, ref = "Normal")
mrna_dds <- DESeq(mrna_dds)
vsd <- varianceStabilizingTransformation(mrna_dds, blind=FALSE)

# Dispersions Plot

plotDispEsts(mrna_dds, main="Dispersion plot")
mrna_res <- results(mrna_dds, name = "Condition_Tumor_vs_Normal")

# MA Plot

plotMA(mrna_res)

# Summing Up

mrna_res_df <- as.data.frame(mrna_res)
mrnaTable <- mrna_res_df
mrnaTable$Gene_id <- rownames(mrnaTable)
summary(mrna_res)
mrna_upreg <- get_upregulated(mrna_res)
mrna_downreg <- get_downregulated(mrna_res)
mrna_counts <- counts(mrna_dds, normalized = T)
mrna_upreg$Gene_id <- rownames(mrna_upreg)
mrna_downreg$Gene_id <- rownames(mrna_downreg)
mrna_res_df$Gene_id <- rownames(mrna_res_df)

#txt Format

write.table(mrna_counts, "C:/Users/thero/Documents/content/Results/mRNAnormcounts.txt", quote = F, sep = "\t")
write.table(mrna_res_df, "C:/Users/thero/Documents/content/Results/mrnaresdeseq2.txt", quote = F, sep = "\t")
write.table(mrna_upreg, "C:/Users/thero/Documents/content/Results/mRNAupreg.txt", quote = F, sep = "\t", row.names = F)
write.table(mrna_downreg, "C:/Users/thero/Documents/content/Results/mRNAdownreg.txt", quote = F, sep = "\t", row.names = F)

#CSV Format

write.csv(as.data.frame(mrna_upreg), file="mrna_upreg.csv")
write.csv(as.data.frame(mrna_downreg), file="mrna_downreg.csv")
write.csv(as.data.frame(mrna_counts), file="mrna_counts.csv")
write.csv(as.data.frame(mrna_res_df), file="mrna_res_df.csv")


# Volcano Plot

keyvals <- ifelse(
    mrna_res_df$log2FoldChange < -1, 'royalblue',
      ifelse(mrna_res_df$log2FoldChange > 1, 'green',
       'grey' ))
  names(keyvals)[keyvals == 'green'] <- 'Up-regulated                        '
  names(keyvals)[keyvals == 'grey'] <- 'Non-significant'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'
EnhancedVolcano(mrna_res_df,
    lab = rownames(mrna_res_df),
    x = 'log2FoldChange',
    y = 'pvalue',
    colCustom = keyvals,
  xlim = c(-5,5),
ylim = c(0, 9), 
 pCutoff = 10e-6,
    pointSize = 2.0,
    labSize = 6.0,
      )
