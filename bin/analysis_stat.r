#R:v4.2.1
#installation des packages
install.packages("BiocManager")
BiocManager::install("DESeq2", force =T )
#BiocManager::install("EnrichmentBrowser")
#BiocManager::install("org.Staphylococcus_aureus.eg.db")
#BiocManager::install("clusterProfiler")
#activer les dépendances
library("DESeq2")
library(ggplot2)
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

# Fichier de nom de gènes : comptage du nombre de ligne du fichier pour ne pas prendre la dernière ligne
file=file("GeneSpecificInformation_COL.tsv", "r")
nbline = 0
while (length(readLines(file, n = 1)) > 0) {
  nbline = nbline + 1
}
close(file)
geneName = read.table("GeneSpecificInformation_NCTC8325.tsv", header=T,sep = "\t",nrows=nbline-1) # on enlève la dernière ligne qui bloque le read.table



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
#plot(log(res$baseMean),res$log2FoldChange,pch=16, col=(res$padj < 0.05)+1 ) # si pval inf alpha : col = 2 (rouge), sinon noir

resdf=data.frame(res)
ggplot(data = resdf, aes(x = log(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), shape = 16,cex=1) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    x = "log(baseMean)",
    y = "log2FoldChange",
    color = "padj < 0.05"
  ) 


##over-representation analysis
gene_list <- row.names(de_genes)
gene_list <- sub("^gene-", "", gene_list)

name_list=geneName[ match(gene_list, geneName[,1]) , 2]

organism_code <- "sao"
kegg_enrich <- enrichKEGG(
  gene = gene_list,
  organism = organism_code,
  keyType = "kegg",
  pvalueCutoff = 0.05
)

#significant_kegg_enrich <- subset(kegg_enrich, padj < 0.05)
#significant_kegg_enrich

