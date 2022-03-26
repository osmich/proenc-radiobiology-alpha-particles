### Setup

# install.packages("tidyverse")
library(tidyverse)
# BiocManager::install("DESeq2", force = TRUE) 
library(DESeq2) 
# BiocManager::install("ashr", force=TRUE)
library(ashr)
library(ggrepel)
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
# BiocManager::install('DOSE')
library(DOSE)

# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library('org.Hs.eg.db')
library(enrichplot)

### CHOIX DE LA METHODE ET DES CONDITIONS DE TEST

 chosen_method = "RNASeq"
 print(paste0("On étudie des données de ",chosen_method))
 chosen_cancer <- "B16" #4T1,B16,MC38 ou MOPC315
 chosen_duration <- "72H" #4H,24H ou 72H
 print(paste0("Etude menée pour le cancer ",chosen_cancer, "et pour une durée d'exposition de ", chosen_duration))


### Load data

if (chosen_method == "RNASeq") {
  data <- read.table("data/RNASeq.dat", header = TRUE, sep="\t",stringsAsFactors = TRUE)
  data = cbind(as.data.frame(data$X), data$X1_MOPC315.BM_ctrl, data$X2_MOPC315.BM_ctrl, data$X6_MOPC315.BM_ctrl, data$X7_MOPC315.BM_ctrl, data$X3_MOPC315.BM_4h.5uCi, data$X8_MOPC315.BM_4h.5uCi, data$X4_MOPC315.BM_24h.5uCi, data$X9_MOPC315.BM_24h.5uCi, data$X5_MOPC315.BM_72h.5uCi, data$X10_MOPC315.BM_72h.5uCi, data$X11_4T1_ctrl, data$X12_4T1_ctrl, data$X16_4T1_ctrl, data$X17_4T1_ctrl, data$X13_4T1_4h.20uCi, data$X18_4T1_4h.20uCi, data$X14_4T1_24h.20uCi, data$X19_4T1_24h.20uCi, data$X15_4T1_72h.20uCi, data$X20_4T1_72h.20uCi, data$X21_B16.F10_ctrl, data$X22_B16.F10_ctrl, data$X26_B16.F10_ctrl, data$X27_B16.F10_ctrl, data$X23_B16.F10_4h.40uCi, data$X28_B16.F10_4h.40uCi, data$X24_B16.F10_24h.40uCi, data$X29_B16.F10_24h.40uCi, data$X25_B16.F10_72h.40uCi, data$X30_B16.F10_72h.40uCi, data$X31_MC.38_ctrl, data$X32_MC.38_ctrl, data$X36_MC.38_ctrl, data$X37_MC.38_ctrl, data$X33_MC.38_4h.40uCi, data$X38_MC.38_4h.40uCi, data$X34_MC.38_24h.40uCi, data$X39_MC.38_24h.40uCi, data$X35_MC.38_72h.40uCi, data$X40_MC.38_72h.40uCi)
  colnames(data) = c("X","X1_MOPC315.BM_ctrl", "X2_MOPC315.BM_ctrl", "X6_MOPC315.BM_ctrl", "X7_MOPC315.BM_ctrl", "X3_MOPC315.BM_4h.5uCi", "X8_MOPC315.BM_4h.5uCi", "X4_MOPC315.BM_24h.5uCi", "X9_MOPC315.BM_24h.5uCi", "X5_MOPC315.BM_72h.5uCi", "X10_MOPC315.BM_72h.5uCi", "X11_4T1_ctrl", "X12_4T1_ctrl", "X16_4T1_ctrl", "X17_4T1_ctrl", "X13_4T1_4h.20uCi", "X18_4T1_4h.20uCi", "X14_4T1_24h.20uCi", "X19_4T1_24h.20uCi", "X15_4T1_72h.20uCi", "X20_4T1_72h.20uCi", "X21_B16.F10_ctrl", "X22_B16.F10_ctrl", "X26_B16.F10_ctrl", "X27_B16.F10_ctrl", "X23_B16.F10_4h.40uCi", "X28_B16.F10_4h.40uCi", "X24_B16.F10_24h.40uCi", "X29_B16.F10_24h.40uCi", "X25_B16.F10_72h.40uCi", "X30_B16.F10_72h.40uCi", "X31_MC.38_ctrl", "X32_MC.38_ctrl", "X36_MC.38_ctrl", "X37_MC.38_ctrl", "X33_MC.38_4h.40uCi", "X38_MC.38_4h.40uCi", "X34_MC.38_24h.40uCi", "X39_MC.38_24h.40uCi", "X35_MC.38_72h.40uCi", "X40_MC.38_72h.40uCi")
  meta <- read.table("./data/RNASeq_meta.txt", header = TRUE, sep=",",stringsAsFactors = TRUE)
}
if (chosen_method == "MassSpec"){
  data <- read.table("data/données_MassSpec_bi.csv", header = TRUE, sep=",")
  meta <- read.table("./data/massSpec_meta.txt", header = TRUE, sep=",",stringsAsFactors = TRUE)
}


### Change columns duration to integer

if (chosen_method == "MassSpec")  {
  data[,c(-1,-2)] <- sapply(data[,c(-1,-2)], as.integer) 
}
if (chosen_method == "RNASeq")  {
  data[,-1] <- sapply(data[,-1], as.integer) 
}

### Setup Desq2 objects


if (chosen_method == "MassSpec")  {
  if (chosen_cancer=="4T1"){
    data_cancer = data[,3:22]
    meta_cancer = meta[1:20,]
  }
  if (chosen_cancer=="B16"){
    data_cancer = data[,23:42]
    meta_cancer = meta[21:40,]
  }
  if (chosen_cancer=="MC38"){
    data_cancer = data[,43:62]
    meta_cancer = meta[41:60,]
  }
  if (chosen_cancer=="MOPC315"){
    data_cancer = data[,63:82]
    meta_cancer = meta[61:80,]
  }
  
}
if (chosen_method == "RNASeq")  {
  if (chosen_cancer=="4T1"){
    data_cancer = data[,12:21]
    meta_cancer = meta[11:20,]
  }
  if (chosen_cancer=="B16"){
    data_cancer = data[,22:31]
    meta_cancer = meta[21:30,]
  }
  if (chosen_cancer=="MC38"){
    data_cancer = data[,32:41]
    meta_cancer = meta[31:40,]
  }
  if (chosen_cancer=="MOPC315"){
    data_cancer = data[,2:11]
    meta_cancer = meta[1:10,]
  }
}

dds <- DESeqDataSetFromMatrix(countData = data_cancer, colData = meta_cancer, design = ~ duration )
if (chosen_method == "MassSpec") {
  rownames(dds) <- data[,2] 
}
if (chosen_method == "RNASeq") {
  rownames(dds) <- data[,1] 
}
dds$duration <- relevel(dds$duration,"CTRL") # necessary to compare all samples against control


### Run analysis

dds <- DESeq(dds)

contrasts <- c("duration", "CTRL", chosen_duration)
res_table_unshrunken <- results(dds, contrast=contrasts , alpha = 0.05)
res_table <- lfcShrink(dds, contrast=contrasts, res=res_table_unshrunken, type="ashr")

normalized_counts <- counts(dds, normalized=TRUE)
data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
normalized_counts <- as.data.frame(normalized_counts)


### Order results by padj values

top20_sig_genes <- as.data.frame(res_table) %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  head(n=20) 		#Extract the first 20 genes

#On enlève les . et les ; des noms de gènes car ils posent problème lors de la jointure ensuite 
row.names(top20_sig_genes) <- gsub(".","", row.names(top20_sig_genes), fixed=TRUE)
row.names(top20_sig_genes) <- gsub(";","", row.names(top20_sig_genes), fixed=TRUE)
row.names(normalized_counts) <- gsub(".","", row.names(normalized_counts), fixed=TRUE)
row.names(normalized_counts) <- gsub(";","", row.names(normalized_counts), fixed=TRUE)

### Normalized counts for top 20 significant genes

top20_sigOE_norm <- subset(normalized_counts, row.names(normalized_counts) %in% row.names(top20_sig_genes))


### save top20 significant genes

write.table(rownames(top20_sig_genes), sprintf("results/top20_genes_%s_%s_%s.txt", chosen_method, chosen_cancer, chosen_duration), row.names=FALSE, col.names=FALSE, quote = FALSE)



############ Volcano plot ##############

## Create a column to indicate which genes to label 

res_tableOE_tb <- as.data.frame(res_table) %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

res_tableOE_tb$genelabels <- row.names(res_tableOE_tb)

ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle(paste("4T1_overexpression",chosen_duration, sep="_")) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 



#### Enhanced Volcano


EnhancedVolcano(toptable = res_table, 
                lab = rownames(res_table), 
                x = 'log2FoldChange', 
                y = 'pvalue',
                title =  paste(paste(chosen_cancer, "CTRL_versus", sep="_"), chosen_duration, sep="_"),
                pCutoff = 10^-15, 
                FCcutoff = 0.8, 
                pointSize = 3.0, 
                labSize = 6.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)


##### Pathway Enrichment analysis

short_res <- res_table[toupper(row.names(res_table)) %in% toupper(rownames(top20_sig_genes)),]
short_res$ENTREZID <- mapIds(org.Hs.eg.db, toupper(row.names(short_res)), 'ENTREZID', 'SYMBOL')

# Biological component

ego <- enrichGO(gene          = short_res$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP", #BP (bioligical component) #CC cellular component
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, # 0.01
                qvalueCutoff  = 0.5, # 0.05
        readable      = TRUE)
head(ego)

goplot(ego)
barplot(ego, showCategory=20)
dotplot(ego, showCategory=30)


# Disease ontology

edo <- enrichDGN(as.vector(short_res$ENTREZID))
head(edo)
barplot(edo, showCategory=20)

dotplot(edo, showCategory=30)
