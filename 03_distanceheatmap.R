



dat <- as.data.frame(assay(rld))
dat <- counts(dds,normalized = T)

sample_order <- c(sample_table$X[c(sample_table$Gender == "F" & sample_table$Condition == "wt,Cre(-)")],
                  sample_table$X[c(sample_table$Gender == "F" & sample_table$Condition == "wt,Cre(+)")],
                  sample_table$X[c(sample_table$Gender == "F" & sample_table$Condition == "Egr1,Cre(-)")],
                  sample_table$X[c(sample_table$Gender == "F" & sample_table$Condition == "Egr1,Cre(+)")],
                  sample_table$X[c(sample_table$Gender == "M" & sample_table$Condition == "wt,Cre(-)")],
                  sample_table$X[c(sample_table$Gender == "M" & sample_table$Condition == "wt,Cre(+)")],
                  sample_table$X[c(sample_table$Gender == "M" & sample_table$Condition == "Egr1,Cre(-)")],
                  sample_table$X[c(sample_table$Gender == "M" & sample_table$Condition == "Egr1,Cre(+)")])

dat <- dat[,sample_order]
sample_table <- sample_table[sample_order,]
colnames(dat) <- paste(sample_table$Condition,sample_table$Gender,sep = ",")

dat <- as.data.frame(dat)
colnames(dat) <- paste(sample_table$sample_id,sample_table$Condition,sample_table$Gender,sep = ",")

dat$Gene <- trans[rownames(dat),]$GeneSymbol
dat <- dat[-which(colnames(dat)=="Gene")]

sampleDists <- dist(t(dat))   #dist默认计算矩阵行与行的距离， 因此需要转置
sampleDistMatrix <- as.matrix(sampleDists)  


column_ha  <- HeatmapAnnotation(Gender = sample_table$Gender,
                                col = list(Gender = c("F" = "red", "M" = "blue")))
row_ha <-  rowAnnotation(Condtion = sample_table$Condition,
                         col = list(Condtion = c("wt,Cre(-)" = "#ffff08","wt,Cre(+)"= "#ff0000",
                                                 "Egr1,Cre(-)" = "#70ad47", "Egr1,Cre(+)" = "#5b9bd5")))
Heatmap(sampleDistMatrix,
        name = "Distance",
        column_title = "Pituitary Gland Distance Heatmap",
        col = c("red","white","blue"),
        top_annotation = column_ha,
        left_annotation = row_ha,
        show_column_names = F,
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(colnames(sampleDistMatrix), rot = 315, just = "left"),
          annotation_height = max_text_width(colnames(sampleDistMatrix))
        ))


tiff(filename = './PIT.Response_to_reviewer/pit.distance.tiff', res = 500, compression = "none",
     width = 9.5, height = 8, units = 'in')

Heatmap(sampleDistMatrix,
        name = "Distance",
        column_title = "Pituitary Gland Distance Heatmap",
        col = c("red","white","blue"),
        cluster_columns = F,
        cluster_rows = F,
        top_annotation = column_ha,
        left_annotation = row_ha,
        show_column_names = F,
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(colnames(sampleDistMatrix), rot = 315, just = "left"),
          annotation_height = max_text_width(colnames(sampleDistMatrix))
        ))

dev.off()