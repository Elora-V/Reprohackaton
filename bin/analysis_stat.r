#R:v4.2.1
#installation des packages
install.packages("BiocManager")
BiocManager::install("DESeq2", force =T )
#BiocManager::install("EnrichmentBrowser")
#BiocManager::install("org.Staphylococcus_aureus.eg.db")
#BiocManager::install("clusterProfiler")
#activer les dépendances
library("DESeq2")
#library("EnrichmentBrowser")
library("KEGGREST")
#library("clusterProfiler")

#analyse statistique de la sortie featureCount
countData = read.table("final_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")

#countData <- read.table("/Users/fionahak/Documents/m2_ami2b/reprohackathon/reprohackathon/results/COUNTING/final_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
#affecter les conditions
colDataFrame <- data.frame(
  Sample = colnames(countData),
  Condition = c("Treatment", "Treatment", "Treatment", "Control", "Control", "Control")
)

##normalization and dispersion estimation
#créer l'objet pour des à partir de la matrice
dds <- DESeqDataSetFromMatrix(countData, colData = colDataFrame, design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds)
#adjust p-values using the Benjamini and Hochberg procedure
res$padj <- p.adjust(res$pvalue, method = "BH")
#genes with an adjusted p value lower than 0.05 were considered differentially expressed
# differentially expressed genes (DEG)
de_genes <- subset(res, padj < 0.05)
#length(de_genes$baseMean) #1487 DEG => 10 de plus que dans l'article

#figure supplémentaire 3
plot(log(res$baseMean),res$log2FoldChange,pch=16, col=(res$padj < 0.05)+1 ) # si pval inf alpha : col = 2 (rouge), sinon noir

##over-representation analysis
gene_list <- row.names(de_genes)
gene_list <- sub("^gene-", "", gene_list)
organism_code <- "sao"
kegg_enrich <- enrichKEGG(
  gene = gene_list,
  organism = organism_code,
  keyType = "kegg",
  pvalueCutoff = 0.05
)

#significant_kegg_enrich <- subset(kegg_enrich, padj < 0.05)
#significant_kegg_enrich

