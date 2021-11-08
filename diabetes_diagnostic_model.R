GSE164416_limma <- DEGs_limma(GSE164416_exp_PCG_log[, GSE164416_cli1$Samples], 
                              GSE164416_cli1[, c("Samples", "Type")],
                              'T2D', 'Healthy')
table(GSE164416_limma$P.Value < 0.05)
GSE164416_limma_filtered <- GSE164416_limma[GSE164416_limma$P.Value < 0.01 & abs(GSE164416_limma$logFC) > log(1.5), ]
dim(GSE164416_limma_filtered)
table(GSE164416_limma_filtered$logFC > 0)
write.csv(GSE164416_limma_filtered,
          file = 'results/GSE164416_limma_filtered.csv')

GSE164416_limma_up <- filter(GSE164416_limma_filtered, logFC > 0) %>% arrange(desc(logFC))
GSE164416_up_genes <- rownames(GSE164416_limma_up)
GSE164416_limma_dn <- filter(GSE164416_limma_filtered, logFC < 0) %>% arrange(logFC)
GSE164416_dn_genes <- rownames(GSE164416_limma_dn)

pdf('PDFs/GSE164416_diffgene_volcano.pdf', width = 5, height = 5)
wb_volcano(rownames(GSE164416_limma), 
           GSE164416_limma$P.Value, 
           GSE164416_limma$logFC, 
           cut_p = 0.01, 
           cut_logFC = log2(1.5),
           xlim = 2, ylab = '-log10(P)')
dev.off()

library(pheatmap)
annotation_col  <- data.frame(Type = factor(GSE164416_cli1$Type))
rownames(annotation_col) <- GSE164416_cli1$Sample
GSE164416_gene_heatmap <- GSE164416_exp[rownames(GSE164416_limma_filtered),
                                        GSE164416_cli1$Sample]
bk=unique(c(seq(-1.2, 1.2, length=100)))
table(GSE164416_cli1$Type)
ann_colors = list(
  Cluster = c(Healthy = '#E64B35', T2D = '#4DBBD5')
)
pdf('PDFs/GSE164416_gene_heatmap.pdf',width = 6,height = 6)
pheatmap(GSE164416_gene_heatmap,
         scale = 'row',
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F, show_colnames = F,
         gaps_col = 18, #gaps_row = 50,
         # cellwidth = 0.8, cellheight = 0.35,
         treeheight_row = 25, treeheight_col = 25,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
         ,annotation_colors =ann_colors
)
dev.off()



GSE164416_up_genes_go_kegg <- clusterProfiler_go_kegg(GSE164416_up_genes,
                                                      title = "UP Genes",
                                                      col = "pvalue")
GSE164416_up_genes_go_kegg$go.kegg.plot

GSE164416_up_genes_go_kegg.plot <- GSE164416_up_genes_go_kegg$go.kegg.plot
GSE164416_up_genes_go_kegg.plot
ggsave(plot = GSE164416_up_genes_go_kegg.plot,
       filename = 'PDFs/GSE164416_up_genes_go_kegg.plot.pdf',
       width = 16.5, height = 8)
GSE164416_up_genes_go_kegg.dat <- GSE164416_up_genes_go_kegg$go.kegg.dat
table(GSE164416_up_genes_go_kegg.dat$pvalue < 0.05, GSE164416_up_genes_go_kegg.dat$TYPE)
write.csv(GSE164416_up_genes_go_kegg.dat,
          file = 'results/GSE164416_up_genes_go_kegg.dat.csv')

GSE164416_dn_genes_go_kegg <- clusterProfiler_go_kegg(GSE164416_dn_genes,
                                                      title = "DOWN Genes",
                                                      col = "pvalue")
GSE164416_dn_genes_go_kegg$go.kegg.plot

GSE164416_dn_genes_go_kegg.plot <- GSE164416_dn_genes_go_kegg$go.kegg.plot
GSE164416_dn_genes_go_kegg.plot
ggsave(plot = GSE164416_dn_genes_go_kegg.plot,
       filename = 'PDFs/GSE164416_dn_genes_go_kegg.plot.pdf',
       width = 19, height = 8)
GSE164416_dn_genes_go_kegg.dat <- GSE164416_dn_genes_go_kegg$go.kegg.dat
table(GSE164416_dn_genes_go_kegg.dat$pvalue < 0.05, GSE164416_dn_genes_go_kegg.dat$TYPE)
write.csv(GSE164416_dn_genes_go_kegg.dat,
          file = 'results/GSE164416_dn_genes_go_kegg.dat.csv')



library(e1071)
library(affy)
library(pROC)
GSE164416_roc_dat <- cbind(t(GSE164416_exp_PCG_log[hub_genes,  GSE164416_cli1$Samples]),
                           GSE164416_cli1[, c("Samples", "Type")])
str(GSE164416_roc_dat)

set.seed(12345)
GSE164416.tObj <- tune.svm(GSE164416_roc_dat[, hub_genes],
                           factor(GSE164416_roc_dat$Type)
                           ,type="C-classification"
                           ,kernel="radial"
                           ,probability = TRUE
                           , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE164416.BestSvm<-GSE164416.tObj$best.model
GSE164416.pre_label <- predict(GSE164416.BestSvm, 
                               GSE164416_roc_dat[, hub_genes], 
                               probability = F,
                               decision.values=T)
GSE164416.pre_label
write.csv(GSE164416.pre_label,
          file = 'files/svm/GSE164416.svm.csv')
dev.off()
GSE164416_svm_res <- cal_auc(GSE164416.pre_label, factor(GSE164416_roc_dat$Type))
table(GSE164416.pre_label,factor(GSE164416_roc_dat$Type))

dir.create('files/svm/')
GSE164416.pre_dat <- read.csv('files/svm/GSE164416.svm.csv')
GSE164416.pre_dat <- merge(GSE164416_roc_dat, GSE164416.pre_dat,
                           by = 'Samples')
GSE164416.pre_dat$Type1 <- as.numeric(factor(GSE164416.pre_dat$Type))

GSE164416_roc <- roc(Type1 ~ Values, 
                     data = GSE164416.pre_dat, auc=TRUE, ci=TRUE)
GSE164416_roc_plot <- ggroc(GSE164416_roc, size = 1.2, legacy.axes = TRUE)
GSE164416_roc_plot <- GSE164416_roc_plot + 
  theme_bw() + labs(title = 'GSE164416') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(GSE164416_roc$auc, 3), '(', 
                         round(GSE164416_roc$ci[1], 3), '-',  
                         round(GSE164416_roc$ci[3], 3), ')'),
           size = 5)
GSE164416_roc_plot
ggsave(plot = GSE164416_roc_plot,
       filename = 'PDFs/GSE164416_roc_plot.pdf',
       width = 5, height = 5)


GSE161355_roc_dat <- cbind(t(GSE161355_exp[hub_genes,  GSE161355_cli1$Samples]),
                           GSE161355_cli1[, c("Samples", "Type")])
str(GSE161355_roc_dat)

set.seed(12345)
GSE161355.tObj <- tune.svm(GSE161355_roc_dat[, hub_genes],
                           factor(GSE161355_roc_dat$Type)
                           ,type="C-classification"
                           ,kernel="radial"
                           ,probability = TRUE
                           , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE161355.BestSvm<-GSE161355.tObj$best.model
GSE161355.pre_label <- predict(GSE161355.BestSvm, 
                               GSE161355_roc_dat[, hub_genes], 
                               probability = F,
                               decision.values=T)
GSE161355.pre_label
write.csv(GSE161355.pre_label,
          file = 'files/svm/GSE161355.svm.csv')
dev.off()
GSE161355_svm_res <- cal_auc(GSE161355.pre_label, factor(GSE161355_roc_dat$Type))
table(GSE161355.pre_label,factor(GSE161355_roc_dat$Type))

dir.create('files/svm/')
GSE161355.pre_dat <- read.csv('files/svm/GSE161355.svm.csv')
GSE161355.pre_dat <- merge(GSE161355_roc_dat, GSE161355.pre_dat,
                           by = 'Samples')
GSE161355.pre_dat$Type1 <- as.numeric(factor(GSE161355.pre_dat$Type))

GSE161355_roc <- roc(Type1 ~ Values, 
                     data = GSE161355.pre_dat, auc=TRUE, ci=TRUE)
GSE161355_roc_plot <- ggroc(GSE161355_roc, size = 1.2, legacy.axes = TRUE)
GSE161355_roc_plot <- GSE161355_roc_plot + 
  theme_bw() + labs(title = 'GSE161355') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(GSE161355_roc$auc, 3), '(', 
                         round(GSE161355_roc$ci[1], 3), '-',  
                         round(GSE161355_roc$ci[3], 3), ')'),
           size = 5)
GSE161355_roc_plot
ggsave(plot = GSE161355_roc_plot,
       filename = 'PDFs/GSE161355_roc_plot.pdf',
       width = 5, height = 5)


GSE156993_roc_dat <- cbind(t(GSE156993_exp[hub_genes,  GSE156993_cli1$Samples]),
                           GSE156993_cli1[, c("Samples", "Type")])
str(GSE156993_roc_dat)

set.seed(12345)
GSE156993.tObj <- tune.svm(GSE156993_roc_dat[, hub_genes],
                           factor(GSE156993_roc_dat$Type)
                           ,type="C-classification"
                           ,kernel="radial"
                           ,probability = TRUE
                           , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE156993.BestSvm<-GSE156993.tObj$best.model
GSE156993.pre_label <- predict(GSE156993.BestSvm, 
                               GSE156993_roc_dat[, hub_genes], 
                               probability = F,
                               decision.values=T)
GSE156993.pre_label
write.csv(GSE156993.pre_label,
          file = 'files/svm/GSE156993.svm.csv')
dev.off()
GSE156993_svm_res <- cal_auc(GSE156993.pre_label, factor(GSE156993_roc_dat$Type))
table(GSE156993.pre_label,factor(GSE156993_roc_dat$Type))

dir.create('files/svm/')
GSE156993.pre_dat <- read.csv('files/svm/GSE156993.svm.csv')
GSE156993.pre_dat <- merge(GSE156993_roc_dat, GSE156993.pre_dat,
                           by = 'Samples')
GSE156993.pre_dat$Type1 <- as.numeric(factor(GSE156993.pre_dat$Type))

GSE156993_roc <- roc(Type1 ~ Values, 
                     data = GSE156993.pre_dat, auc=TRUE, ci=TRUE)
GSE156993_roc_plot <- ggroc(GSE156993_roc, size = 1.2, legacy.axes = TRUE)
GSE156993_roc_plot <- GSE156993_roc_plot + 
  theme_bw() + labs(title = 'GSE156993') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(GSE156993_roc$auc, 3), '(', 
                         round(GSE156993_roc$ci[1], 3), '-',  
                         round(GSE156993_roc$ci[3], 3), ')'),
           size = 5)
GSE156993_roc_plot
ggsave(plot = GSE156993_roc_plot,
       filename = 'PDFs/GSE156993_roc_plot.pdf',
       width = 5, height = 5)


GSE163980_roc_dat <- cbind(t(GSE163980_exp[hub_genes,  GSE163980_cli1$Samples]),
                           GSE163980_cli1[, c("Samples", "Type")])
str(GSE163980_roc_dat)

set.seed(12345)
GSE163980.tObj <- tune.svm(GSE163980_roc_dat[, hub_genes],
                           factor(GSE163980_roc_dat$Type)
                           ,type="C-classification"
                           ,kernel="radial"
                           ,probability = TRUE
                           , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE163980.BestSvm<-GSE163980.tObj$best.model
GSE163980.pre_label <- predict(GSE163980.BestSvm, 
                               GSE163980_roc_dat[, hub_genes], 
                               probability = F,
                               decision.values=T)
GSE163980.pre_label
write.csv(GSE163980.pre_label,
          file = 'files/svm/GSE163980.svm.csv')
dev.off()
GSE163980_svm_res <- cal_auc(GSE163980.pre_label, factor(GSE163980_roc_dat$Type))
table(GSE163980.pre_label,factor(GSE163980_roc_dat$Type))

dir.create('files/svm/')
GSE163980.pre_dat <- read.csv('files/svm/GSE163980.svm.csv')
GSE163980.pre_dat <- merge(GSE163980_roc_dat, GSE163980.pre_dat,
                           by = 'Samples')
GSE163980.pre_dat$Type1 <- as.numeric(factor(GSE163980.pre_dat$Type))

GSE163980_roc <- roc(Type1 ~ Values, 
                     data = GSE163980.pre_dat, auc=TRUE, ci=TRUE)
GSE163980_roc_plot <- ggroc(GSE163980_roc, size = 1.2, legacy.axes = TRUE)
GSE163980_roc_plot <- GSE163980_roc_plot + 
  theme_bw() + labs(title = 'GSE163980') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(GSE163980_roc$auc, 3), '(', 
                         round(GSE163980_roc$ci[1], 3), '-',  
                         round(GSE163980_roc$ci[3], 3), ')'),
           size = 5)
GSE163980_roc_plot
ggsave(plot = GSE163980_roc_plot,
       filename = 'PDFs/GSE163980_roc_plot.pdf',
       width = 5, height = 5)
