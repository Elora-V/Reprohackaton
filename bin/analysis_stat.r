args <- commandArgs(trailingOnly = TRUE)

#arguments
final_count_matrix <- args[1]
GSE139659_IPvsctrl.complete <- args[2]
GeneSpecificInformation_NCTC8325 <- args[3]
geneTranslation <- args[4]

#data = 0 si nos données, 1 sinon

for (data in c(1,0)){

  #R:v4.2.1
  #activer les dépendances
  library("DESeq2")
  library(ggplot2)
  library(ggrepel)


  ### Counting matrix
  if (data == 0){
    #analyse statistique de la sortie featureCount
    countData = read.table(final_count_matrix, header = TRUE, row.names = 1, sep = "\t")
    #affecter les conditions
    colDataFrame <- data.frame(
      Sample = colnames(countData),
      Condition = c("Treatment", "Treatment", "Treatment", "Control", "Control", "Control")
    )
  } else {
    # pr celle de l'article
    countData=read.delim(GSE139659_IPvsctrl.complete,header=T,row.names = NULL) # pour tester avec la vrai matrice
    countData=countData[1:2967,]
    rownames(countData)=countData$Name
    countData=countData[,6:11]
    countData= na.omit(countData)
    colDataFrame <- data.frame(
      Sample = colnames(countData),
      Condition = c("Control", "Control", "Control","Treatment", "Treatment", "Treatment")
    )
  }

  ### Genes names matrix

  # Fichier de nom de gènes : comptage du nombre de ligne du fichier pour ne pas prendre les dernières lignes
  file=file(GeneSpecificInformation_NCTC8325, "r")
  nbline = 0
  while (length(readLines(file, n = 1)) > 0) {
    nbline = nbline + 1
  }
  close(file)
  geneName = read.table(GeneSpecificInformation_NCTC8325, header=T,sep = "\t",nrows=nbline-3) # on enlève les dernières lignes qui bloque le read.table

  ### Differential expression analysis

  ##normalization and dispersion estimation
  #créer l'objet pour deseq2 à partir de la matrice
  dds <- DESeqDataSetFromMatrix(countData, colData = colDataFrame, design = ~ Condition)
  dds <- DESeq(dds)
  res <- results(dds)
  #Changement des noms de lignes et Ajout des symboles des genes du tableau de résultat
  rownames(res)=sub("^gene-", "", rownames(res))
  name=geneName[ match(rownames(res),geneName[,1]) , 2]
  indexNa=which(is.na(name))
  name[indexNa]=rownames(res)[indexNa]
  res$name=name

  #adjust p-values using the Benjamini and Hochberg procedure
  res$padj <- p.adjust(res$pvalue, method = "BH")
  #genes with an adjusted p value lower than 0.05 were considered differentially expressed
  # differentially expressed genes (DEG)
  de_genes <- subset(res, padj < 0.05)
  #length(de_genes$baseMean) #1487 DEG => 10 de plus que dans l'article

  if (data == 0){
    RES=res
    GeneDiff=de_genes
  } else {
    RES_A=res
    GeneDiff_A=de_genes
  }

  #### Figure supplémentaire 3 : MA plot

  res$log2FoldChange=ifelse(res$log2FoldChange > 4, 4, ifelse(res$log2FoldChange < -4, -4, res$log2FoldChange))
  resdf=data.frame(res)

  maplot=ggplot(data = resdf, aes(x = log(baseMean), y = log2FoldChange)) +
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

  transl=read.table(geneTranslation,sep="\t") # separateur ne sépare rien, mais c'est normal car le format des lignes diffèrent
  transl[,1] <- sub("^\\s+", "", transl[,1])  # supprime les premiers espaces
  transl[,1] <- sub("^(\\S+).*", "\\1", transl[,1]) # on récupère juste le nom du gène
  # sub(pattern, replacement, data)
  # ^ indique que c'est en début de chaine
  # (\\S+) correspond à un ou plusieurs caractères consécutifs qui ne sont pas des espaces blancs
  # "\1" fait référence au premier groupe de capture
  transl[,1] <- sub("^sao:", "", transl[,1])

  ## result for translation

  index_transl=match(transl[,1],rownames(res))
  index_transl=na.omit(index_transl) # 3 NA
  res_transl=data.frame(res[index_transl,])

  if (data == 0){
    REST=res_transl
  } else {
    REST_A=res_transl
  }

  # MA plot (translation)

  maplot_transl <- ggplot(data = res_transl, aes(x = log2(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05), shape = 16, cex = 1.5) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "log(baseMean)",
      y = "log2FoldChange",
      color = "padj < 0.05",
      title = "MA-plot of translation genes"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_label_repel(data = subset(res_transl, !grepl("^SAOUHSC", name)),
                     aes(label = name),
                     box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.3, "lines"))

  print(maplot_transl)


  if (data == 0){
    titre="MA-plot.pdf"
  } else {
    titre="MA-plot_article.pdf"
  }
  pdf(titre, width = 10, height = 7)
  print(maplot)            # in the first page of PDF
  print(maplot_transl)     # in the second page of the PDF
  dev.off()
}

