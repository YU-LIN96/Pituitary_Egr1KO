rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup=c("Condition", "Gender"),ntop = 3000, returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))

pca <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Gender)) + 
  geom_point(size=3) + 
  scale_shape_manual(values=c(16,15)) +
  scale_color_manual(values=c('#ffff08','#ff0000', '#70ad47','#5b9bd5')) +
  ggtitle("Pituitary Gland PCA") + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2,fill = NA))
pca
# size = 600*500
ggsave('./PIT.Response_to_reviewer/pit.pca.tiff', device='tiff', dpi=500, 
       height = 5, width = 6, unit = 'in')