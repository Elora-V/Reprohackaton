setwd("~/Documents/M2_AMI2B/Reprohackaton/reprohackathon/bin")
rm(list = ls())

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
library("EnrichmentBrowser")
library("KEGGREST")
library("clusterProfiler")


### Counting matrix 

#analyse statistique de la sortie featureCount
countData = read.table("final_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")

#countData <- read.table("/Users/fionahak/Documents/m2_ami2b/reprohackathon/reprohackathon/results/COUNTING/final_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
#affecter les conditions
colDataFrame <- data.frame(
  Sample = colnames(countData),
  Condition = c("Treatment", "Treatment", "Treatment", "Control", "Control", "Control")
)

### Genes names matrix

# Fichier de nom de gènes : comptage du nombre de ligne du fichier pour ne pas prendre les dernières lignes
file=file("GeneSpecificInformation_NCTC8325.tsv", "r")
nbline = 0
while (length(readLines(file, n = 1)) > 0) {
  nbline = nbline + 1
}
close(file)
geneName = read.table("GeneSpecificInformation_NCTC8325.tsv", header=T,sep = "\t",nrows=nbline-3) # on enlève les dernières lignes qui bloque le read.table


### Differential expression analysis

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

# get name of gene

name=geneName[ match(rownames(res),geneName[,1]) , 2]
indexNa=which(is.na(name))
name[indexNa]=rownames(res)[indexNa]
res$name=name

#### Figure supplémentaire 3 : MA plot

res$log2FoldChange=ifelse(res$log2FoldChange > 4, 4, ifelse(res$log2FoldChange < -4, -4, res$log2FoldChange))
resdf=data.frame(res)

ggplot(data = resdf, aes(x = log(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), shape = 16,cex=1) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    x = "log(baseMean)",
    y = "log2FoldChange",
    color = "padj < 0.05"
  )  +
  scale_y_continuous(limits = c(-4, 4)) +
geom_hline(yintercept = 0, linetype = "dashed") +
   labs(title = "MA-plot of complete RNA-seq dataset")+
  theme(plot.title = element_text(hjust = 0.5))


#### Figure 3c : MA plot for translation

# genes for translation

transl=read.table("geneTranslation.txt",sep="\t") # separateur ne sépare rien, mais c'est normal car le format des lignes diffèrent
transl[,1] <- sub("^\\s+", "", transl[,1])  # supprime les premiers espaces
transl[,1] <- sub("^(\\S+).*", "\\1", transl[,1]) # on récupère juste le nom du gène
# sub(pattern, replacement, data)
# ^ indique que c'est en début de chaine
# (\\S+) correspond à un ou plusieurs caractères consécutifs qui ne sont pas des espaces blancs
# "\1" fait référence au premier groupe de capture
transl[,1] <- sub("^sao:", "", transl[,1])

# verification avec un autre fichier qu'on a les bons gènes (à supprimer plus tard)
Transl2=read.table("geneTranslation2.txt",sep="\t")
Transl2[,1] <- sub("^\\s+", "", Transl2[,1])  # supprime les premiers espaces
Transl2[,1] <- sub("^(\\S+).*", "\\1", Transl2[,1]) # recupère le premeir mot

## result for translation

rownames(res)=sub("^gene-", "", rownames(res))
index_transl=match(transl[,1],rownames(res))
index_transl=na.omit(index_transl) # 3 NA
res_transl=data.frame(res[index_transl,])



# MA plot (traslation)


ggplot(data = res_transl, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), shape = 16,cex=2) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    x = "log(baseMean)",
    y = "log2FoldChange",
    color = "padj < 0.05"
  )  +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "MA-plot of translation genes")+
  theme(plot.title = element_text(hjust = 0.5))



# 
# organism_code <- "sao"
# kegg_enrich <- enrichKEGG(
#   gene = gene_list,
#   organism = organism_code,
#   keyType = "kegg",
#   pvalueCutoff = 0.05
# )

#significant_kegg_enrich <- subset(kegg_enrich, padj < 0.05)
#significant_kegg_enrich

