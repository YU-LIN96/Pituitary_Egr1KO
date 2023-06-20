res_1 <- results(dds,
                 contrast = c('condition',"ffPos_F","ffPos_M"))
res_1 <- as.data.frame(res_1)
res_1 <- res_1[!is.na(res_1$padj),]
res_1$gene_symbol <- trans[rownames(res_1),]$GeneSymbol
res_1 <- res_1[!is.na(res_1$gene_symbol),]
# write.csv(res_1, file = "./PIT_DE_Results/Egr1(ff)Cre(+)F_vs_M.csv")

res_1$padj[res_1$padj < 10e-30] <- 10e-30

keyvals <- ifelse(
  res_1$log2FoldChange < -log2(1.5) & res_1$padj < 0.05, 'blue',
  ifelse(res_1$log2FoldChange > log2(1.5) & res_1$padj < 0.05, 'red',
         'black'))
names(keyvals)[keyvals == 'red'] <- 'Upregulated'
names(keyvals)[keyvals == 'black'] <- 'N.S.'
names(keyvals)[keyvals == 'blue'] <- 'Downregulated'


p1 <- EnhancedVolcano(res_1,
                      lab = res_1$gene_symbol,
                      labSize = 8,
                      selectLab = res_1$gene_symbol[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                      colCustom = keyvals,
                      x = 'log2FoldChange',
                      y = 'padj',
                      # ylab = bquote(~Log[10]~ 'Adjusted p'),
                      caption = "",
                      ylab = "",
                      xlab = "",
                      pCutoff = 0.05,
                      FCcutoff = log2(1.5),
                      subtitle = "",
                      # title = bquote("Egr1"^"ff"~"Nes"^"Cre+"),
                      title = "",
                      legendPosition = 'none',
                      colAlpha = 1,
                      gridlines.major = T,
                      gridlines.minor = F) +
  ggplot2::coord_cartesian(xlim=c(-15, 15),ylim = c(0,30)) +
  ggplot2::scale_x_continuous(
    breaks=c(-15,-10,-5,-1,1,5,10,15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size = 2,fill = NA))

ggsave('./PIT.Response_to_reviewer/pit.p1.tiff', device='tiff', dpi=500, 
       height = 8, width = 10, unit = 'in')