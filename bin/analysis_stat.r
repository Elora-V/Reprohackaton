args <- commandArgs(trailingOnly = TRUE)

#arguments
final_count_matrix <- args[1]
GSE139659_IPvsctrl.complete <- args[2]
GeneSpecificInformation_NCTC8325 <- args[3]
geneTranslation <- args[4]

#R:v4.2.1
#dependencies activation
library("DESeq2")
library(ggplot2)
library(ggrepel)


# We execute the scripts for the data of the article (data=1), then our data (data=0)
for (data in c(1,0)){

  ### Opening Counting matrix ----
  if (data == 0){ # our data
    # output of featureCount
    countData = read.table(final_count_matrix, header = TRUE, row.names = 1, sep = "\t")
    colnames(countData)=c("Treatment1","Treatment2","Treatment3","Control1","Control2","Control3")
    rownames(countData)=sub("^gene-", "", rownames(countData))
    DATA=countData
    # association of condition to colums (replicates)
    colDataFrame <- data.frame(
      Sample = colnames(countData),
      Condition = c("Treatment", "Treatment", "Treatment", "Control", "Control", "Control")
    )
  } else {
    # data of article
    tab=read.delim(GSE139659_IPvsctrl.complete,header=T,row.names = NULL) # pour tester avec la vrai matrice
    tab=tab[1:2967,]
    rownames(tab)=tab$Name
    Deseq2_A=tab[,12:17]
    countData=tab[,6:11]
    countData= na.omit(countData)
    DATA_A=countData
    colDataFrame <- data.frame(
      Sample = colnames(countData),
      Condition = c("Control", "Control", "Control","Treatment", "Treatment", "Treatment")
    ) 
    
  }

  ### Opening Genes names matrix ----

  # Gene name file : counting of the number of lines of the file, needed to not take into account the last lines (error)
  file=file(GeneSpecificInformation_NCTC8325, "r")
  nbline = 0
  while (length(readLines(file, n = 1)) > 0) {
    nbline = nbline + 1
  }
  close(file)
  # opening the file as the table
  geneName = read.table(GeneSpecificInformation_NCTC8325, header=T,sep = "\t",nrows=nbline-3) 

  
  
  # Creation of deseq2 objet from the matrix
  dds <- DESeqDataSetFromMatrix(countData, colData = colDataFrame, design = ~ Condition)
  
  
  ### Exploratory analysis ----
  
  # Histogram of read counts
  histCounts=ggplot( log2(data.frame(x=c(counts(dds)))), aes(x = x)) +
    geom_histogram() +
    labs(title = "Histogram of read counts (log2)",
         x = "Log2(Read Counts)",
         y = "Frequency")
  
  
  # Number of read per replicates
  countlog=log2(counts(dds)+1)
  n=dim(countlog)[1]
  if (data == 0){
    legend=c(rep("Treatment",3*n), rep("Control",3*n))
    label=c(rep("Treatment1",n),rep("Treatment2",n),rep("Treatment3",n),rep("Control1",n),rep("Control2",n),rep("Control3",n))
    
  }else{
    legend=c(rep("Control",3*n), rep("Treatment",3*n))
    label=c(rep("Control1",n),rep("Control2",n),rep("Control3",n),rep("Treatment1",n),rep("Treatment2",n),rep("Treatment3",n))
  }
  df=data.frame(x=c(countlog),  y=label, col=legend )
  boxplotCounts=ggplot(df, aes(x=x, y=y, col = legend)) +
    geom_boxplot()
  
  
  # Correlation heatmap between replicates
  cor= cor(counts(dds))
  # heatmap in the pdf directly  
  
  # PCA
  # VST transformation (variance stabilizing transformation)
  vsdata = vst(dds, blind=FALSE)
  
  # get PCA coordinates 
  PCA = plotPCA(vsdata, intgroup = "Condition", returnData = TRUE)
  
  # plot PCA
  PCAplot=ggplot(PCA, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point()+
    geom_label_repel(data = PCA ,  
                     aes(label = rownames(PCA)),
                     box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.3, "lines")) 
  
  ### Differential expression analysis ----
  
  # normalization and dispersion estimation
  dds <- DESeq(dds)
  res <- results(dds)
  
  # get symbol of genes and add them in the deseq2 results
  rownames(res)=sub("^gene-", "", rownames(res))
  name=geneName[ match(rownames(res),geneName[,1]) , 2]
  indexNa=which(is.na(name))
  name[indexNa]=rownames(res)[indexNa] # when there is no match in the gene name table, we use the full name
  res$name=name # we add a column with the short names (or symbols)

  # adjusting the p-values using the Benjamini and Hochberg procedure
  res$padj <- p.adjust(res$pvalue, method = "BH")
  # genes with an adjusted p value lower than 0.05 were considered differentially expressed genes (DEG)
  de_genes <- subset(res, padj < 0.05)
  #length(de_genes$baseMean) #1487 DEG => 10 de plus que dans l'article

  # We keep the results for those data in order to compare the result with the two dataset
  if (data == 0){
    RES=res # result of deseq2
    GeneDiff=de_genes # DEG
  } else {
    RES_A=res # result of deseq2 article
    GeneDiff_A=de_genes # DEG article
  }

  #### Figure 3 : full MA plot ----

  # In the article, the y limit of the plot are -4 and 4 :
  res$log2FoldChange=ifelse(res$log2FoldChange > 4, 4, ifelse(res$log2FoldChange < -4, -4, res$log2FoldChange))
  resdf=data.frame(res) # conversion in dataframe

  # full MA-plot with ggplot2
  maplot=ggplot(data = resdf, aes(x = log(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05), shape = 16,cex=1) + # the color depend on the pvalue
    scale_color_manual(values = c("black", "red")) +
    labs(
      x = "log(baseMean)",
      y = "log2FoldChange",
      color = "padj < 0.05"
    )  +
    scale_y_continuous(limits = c(-4, 4)) +
    geom_hline(yintercept = 0, linetype = "dashed") + # line at 0
    labs(title = "MA-plot of complete RNA-seq dataset")+
    theme(plot.title = element_text(hjust = 0.5))


  #### Figure 3c : MA plot for translation ----

  # get list of genes for translation

  transl=read.table(geneTranslation,sep="\t") # not a good format, all the informations on the same line
  transl[,1] <- sub("^\\s+", "", transl[,1])  # suppress first spaces 
  transl[,1] <- sub("^(\\S+).*", "\\1", transl[,1]) # get gene name (the first word)
  # Explanation of the command :
  # sub(pattern, replacement, data)
  # ^ : beggining of the string
  # (\\S+) : one or more consecutives characters which aren't spaces
  # "\1" : first element
  transl[,1] <- sub("^sao:", "", transl[,1]) # suppress "sao:" of the name

  ## get result for translation

  index_transl=match(transl[,1],rownames(res))
  index_transl=na.omit(index_transl) # 3 NA
  res_transl=data.frame(res[index_transl,]) # we keep the dese2 results only for the translation genes

  if (data == 0){
    REST=res_transl # deseq2 results for translation
  } else {
    REST_A=res_transl  # deseq2 results for translation, for the article
  }

  # MA plot translation

  res_transl["SAOUHSC_00475","name"]="pth" # one gene with wrong symbol
  
  maplot_transl = ggplot(data = res_transl, aes(x = log2(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05), shape = 16, cex = 1.5) + # color depend on pvalue
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "log(baseMean)",
      y = "log2FoldChange",
      color = "padj < 0.05",
      title = "MA-plot of translation genes"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_label_repel(data = res_transl[match(c("frr","infA","infB","infC","pth","tsf"), res_transl$name),] ,  # get only gene in the list given
                     aes(label = name),
                     box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.3, "lines"))
  
  
  #### PDF ----
  # title of the pdf with the MA plot
  if (data == 0){
    titre="MA-plot.pdf"
  } else {
    titre="MA-plot_article.pdf"
  }
  pdf(titre, width = 10, height = 7)
  print(histCounts)
  print(boxplotCounts)
  print(heatmap(cor))
  print(PCAplot)
  print(hist(res$padj))
  print(maplot)     
  print(maplot_transl)    
  dev.off() 
  
}

#### Comparing results ----

print("################")

### Comparing the counting matrix
print(paste("There is ",dim(RES)[1],"genes (our data)."))
#  "There is  2967 genes (our data)."
print(paste("There is ",dim(RES_A)[1],"genes (data article)."))
#  "There is  2967 genes (data article)."

# change order line to be the same as the results in the article
orderGene=match(rownames((RES_A)), rownames(RES))
print(paste("There is ",length(which(is.na(orderGene)==TRUE)), "genes in the article's data not in our data."))
# "There is  0 genes in the article's data not in our data."

# => Same genes in both matrix

# Data with article's data and our data

DATA=DATA[orderGene,]
fullCounts=cbind(DATA,DATA_A)
colDataFrame= data.frame(
  Sample = colnames(fullCounts),
  Condition = c(rep("Treatment",3),rep( "Control",6),rep("Treatment",3))
)
ddsfull= DESeqDataSetFromMatrix(fullCounts, colData = colDataFrame, design = ~ Condition)
vsdata = vst(ddsfull, blind=FALSE)
PCA = plotPCA(vsdata, intgroup = "Condition", returnData = TRUE)
PCAplot=ggplot(PCA, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point()+
  geom_label_repel(data = PCA ,  
                   aes(label = rownames(PCA)),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.3, "lines")) 
cor= as.matrix(cor(counts(ddsfull)))

# => Treatments of both dataset are closer together than to the control, same idea for control.
# => Treatments1 closer to IP1 (treatment 1 for article), same for the others.
## => Both counting matrix are similar.


### Comparing deseq2 results

# DESeq2 article vs DESeq2 conting matrix article

Deseq2_A=Deseq2_A[match(rownames(RES_A),rownames(Deseq2_A)),]
fullRES=data.frame(fold1=Deseq2_A$log2FoldChange,fold2=RES_A$log2FoldChange)
fullRES=na.omit(fullRES)
dim(fullRES) #2736
foldSimilaire=sum(abs(fullRES$fold1-fullRES$fold2)<= 10**(-7) )
print(paste("Il y a ",foldSimilaire," gènes qui ont leurs foldchanges égaux pour une précision de 10**-7, entre les deux formes de résultats de l'article parmi ",dim(fullRES)[1]," gènes non NA"))
# "Il y a  2736  gènes qui ont leurs foldchanges égaux pour une précision de 10**-7,
# entre les deux formes de résultats de l'article parmi  2736  gènes non NA" 

# => ~ Meme logFoldChange entre les deux types de données de l'article


# DEseq between our data and article's data

RES=RES[match(rownames(RES_A),rownames(RES)),]
fullRES=data.frame(fold1=RES$log2FoldChange,fold2=RES_A$log2FoldChange)
fullRES=na.omit(fullRES)
dim(fullRES) #2730
foldSimilaire=sum(abs(fullRES$fold1-fullRES$fold2)<= 10**(-1) )
print(paste("Il y a ",foldSimilaire," gènes qui ont leurs foldchanges différe au maximum de 0.1, entre les résultats de l'article et les nôtres parmi ",dim(fullRES)[1]," gènes non NA"))
# "Il y a  1677  gènes qui ont leurs foldchanges différe au maximum de 0.1,
# entre les résultats de l'article et les nôtres parmi  2730  gènes non NA"

foldSimilaire=sum(abs(fullRES$fold1-fullRES$fold2)<= 0.2)
print(paste("Il y a ",foldSimilaire," gènes qui ont leurs foldchanges différe au maximum de 0.2, entre les résultats de l'article et les nôtres parmi",dim(fullRES)[1]," gènes non NA"))
"Il y a  2116  gènes qui ont leurs foldchanges différe au maximum de 0.2, 
# entre les résultats de l'article et les nôtres parmi 2730  gènes non NA"

moy=mean(abs(fullRES$fold1-fullRES$fold2) )
summary(c(fullRES[,1]),fullRES[,2])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-5.04344 -0.94906 -0.04698  0.04562  0.92211  5.88998 
print(paste("En moyenne diffère de ",round(moy,2)," en foldchange pour des valeurs de foldchange entre -5 et 5, avec pour moyenne 0.046."))

# =>  "En moyenne diffère de  0.25  en foldchange pour des valeurs de foldchange entre -5 et 5, avec pour moyenne 0.046."

# DEG

GeneDiff_deseq= subset(Deseq2_A, padj < 0.05)

print(paste("There is ",dim(GeneDiff)[1]," DEG genes (our data)."))
#"There is  1487  DEG genes (our data)."
print(paste("There is ",dim(GeneDiff_A)[1],"DEG genes (data article)."))
# "There is  1479 DEG genes (data article)."
print(paste("There is ",dim(GeneDiff_deseq)[1],"DEG genes (data article deseq2)."))
# "There is  1479 DEG genes (data article deseq2)."


CommunDEG=length(intersect(rownames((GeneDiff_A)), rownames(GeneDiff)))
CommunDEG_A=length(intersect(rownames((GeneDiff_A)), rownames(GeneDiff_deseq)))

print(paste("Il y a ",CommunDEG_A," genes différentiels commun entre les deux jeux de l'article."))
# "Il y a  1479  genes différentiels commun entre les deux jeux de l'article."
print(paste("Il y a ",CommunDEG," genes différentiels commun entre l'article et nos résultats."))
# "Il y a  1352  genes différentiels commun entre l'article et nos résultats."

### Comparing differential translation gene

print(paste("Il y a ",dim(REST)[1], "gene de traduction parmi les gènes étudiés."))
# "Il y a  167 gene de traduction parmi les gènes étudiés."
DEGtrans=subset(REST,padj<0.05)
DEGtrans_A=subset(REST_A,padj<0.05)
print(paste("Il y a ", dim(DEGtrans)[1], " genes de traduction exprimés différentiellement dans nos données."))
#"Il y a  57  genes de traduction exprimés différentiellement dans nos données."
print(paste("Il y a ", dim(DEGtrans_A)[1], " genes de traduction exprimés différentiellement dans nos données."))
# Il y a  53  genes de traduction exprimés différentiellement dans nos données."

CommunDEGTrans=length(intersect(rownames(DEGtrans),rownames(DEGtrans_A)))
print(paste("Il y a ",CommunDEGTrans," genes  de traductions communs exprimés différentiellement entre nos données et l'article."))
#  "Il y a  47  genes de traductions communs exprimés différentiellement entre nos données et l'article."



# PDF ----
pdf("Comparaison_resultats.pdf", width = 10, height = 7)
print(PCAplot)
print(heatmap(cor))
dev.off()



