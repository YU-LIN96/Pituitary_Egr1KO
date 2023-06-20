```{r}
library(glmmSeq)
library(emmeans)
library(volcano3D)
library(plotly)
library(fmsb)
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
sample_table$phenotype <- ifelse(
  sample_table$Condition == "Egr1,Cre(+)", 'KO',
  'Normal')

a <- data.frame(matrix(unlist(str_split(sample_table$Condition, ",")), ncol = 2, byrow=TRUE))

sample_table$genotype <- a$X1
sample_table$Cre <- ifelse(
  a$X2 == "Cre(+)", 'Pos',
  'Neg')

# sample_table <- sample_table[c(3,5:9)]
```


```{r}
trans <- read.delim("X:/Xie Lab/Data/FA_bulk/rna_featurecount/Mus_musculus.GRCm39.107_gene_annotation_table.txt")
rownames(trans) <- trans$gene_id
```


```{r}
cont <- raw_counts
cont$gene_symbol <- trans[rownames(cont),]$GeneSymbol
cont <- cont[!duplicated(cont$gene_symbol),]
cont <- cont[!is.na(cont$gene_symbol),]
cont <- cont[!str_detect(cont$gene_symbol, pattern = "gene_source"),]
rownames(cont) <- cont$gene_symbol
cont <- cont[1:24]
```

```{r}
dds <- DESeqDataSetFromMatrix(countData = cont,
                              colData =  sample_table, 
                              design = ~condition)
dds <- DESeq(dds)

disp <- setNames(dispersions(dds), rownames(cont))
sizeFactor <- estimateSizeFactorsForMatrix(cont)
```



```{r}
glmm_res <- glmmSeq(~ phenotype+genotype+Cre+Gender+ (1|sample_id),
                    countdata = cont,
                    metadata = sample_table,
                    dispersion = disp,
                    sizeFactors = sizeFactor,
                    progress = T,
                    cores = 10)
summary(glmm_res, gene = "Egr1")
save(glmm_res, file = "glmm_res.Rdata")

write.csv(summary(glmm_res),file = "mix.model.results.csv")

```