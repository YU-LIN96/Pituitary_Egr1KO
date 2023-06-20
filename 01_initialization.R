```{r}
library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggdendro)
library(circlize)
library(EnhancedVolcano)
library(ggvenn)
library(reshape2)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
```



```{r}
raw_counts <- read.csv("./PIT/c_union_rawCounts_newSamples.csv",
                       header = T, sep = ',', stringsAsFactors = F, row.names = 1)

sample_table <- read.csv("./PIT/c_sampletable.csv")

rownames(sample_table) <- colnames(raw_counts)

detail_inf <- data.frame(matrix(unlist(str_split(sample_table$X, "_")), ncol = 3, byrow=TRUE))
colnames(detail_inf) <- c("sample_id","Exp_condition","Gender")
sample_table <- cbind(sample_table,detail_inf)
sample_table$condition <- paste(sample_table$Exp_condition,sample_table$Gender,sep = "_")
sample_table$Exp_condition <- factor(sample_table$Exp_condition, levels = c("wtNeg","wtPos","ffNeg","ffPos"))


sample_table$Condition[sample_table$Exp_condition == "wtNeg"] <- "wt,Cre(-)"
sample_table$Condition[sample_table$Exp_condition == "wtPos"] <- "wt,Cre(+)"
sample_table$Condition[sample_table$Exp_condition == "ffNeg"] <- "Egr1,Cre(-)"
sample_table$Condition[sample_table$Exp_condition == "ffPos"] <- "Egr1,Cre(+)"

sample_table$Condition <- factor(sample_table$Condition, levels = c("wt,Cre(-)","wt,Cre(+)","Egr1,Cre(-)","Egr1,Cre(+)"))
```



```{r}
trans <- read.delim("X:/Xie Lab/Data/FA_bulk/rna_featurecount/Mus_musculus.GRCm39.107_gene_annotation_table.txt")
rownames(trans) <- trans$gene_id
trans.length <- as.data.frame(t(as.data.frame(strsplit(trans$Chromosome, split = ":"))))
trans.length <- as.data.frame(t(as.data.frame(strsplit(trans.length$V2, split = "-"))))
trans$length <- as.numeric(trans.length$V2)-as.numeric(trans.length$V1)
```


```{r}
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData =  sample_table, 
                              design = ~condition)
dds <- DESeq(dds)
```