#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####  SET WORKING DIRECTORY ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
options(stringsAsFactors = FALSE)
setwd("/Users/k2367592/Desktop/Ligand_Data/results-S3/")
require("DESeq2");require("ggplot2");require("ggrepel")
library("DGEobj.utils")

### Gene labels
require("biomaRt")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_all_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)


####  LOAD DATA ####
#~~~~~~~~~~~~~~~~~~#

counts <- read.table("featurecounts/refined_counts_mat.csv", sep = ",", header = TRUE)
# counts <- read.table("featurecounts/RNA_seq_featurecounts-arch.txt", sep = "\t", header = TRUE)

# Remove columns which we don't use, and sort
# order = sort(names(counts)[!grepl('\\.', names(counts))])
order = sort(names(counts)[!grepl('Length', names(counts))])
geneLength = counts[,c('Geneid', 'Length')]
counts <- counts[,order]

# Rename by substitution
order_rename = gsub('X1', 'MCF10A_', gsub('S', 'SUM159_', gsub('M', 'MDAMB231_', order)))
order_rename = gsub('_C', '_ctrl', gsub('_D', '_dox', order_rename))
colnames(counts) = order_rename

# Correct the rotational mis-labelling issue
colnames(counts)[colnames(counts) == 'MDAMB231_ctrl2'] <- 'SUM159_ctrl2-a'
colnames(counts)[colnames(counts) == 'SUM159_ctrl2'] <- 'MCF10A_ctrl2-a'
colnames(counts)[colnames(counts) == 'MCF10A_ctrl2'] <- 'MDAMB231_ctrl2-a'
colnames(counts) = gsub('-a', '', colnames(counts))

rownames(counts) <- counts$Geneid
rownames(geneLength) <- counts$Geneid
counts <- subset(counts, select = -Geneid)
geneLength <- subset(geneLength, select = -Geneid)
counts <- counts[,sort(colnames(counts))]
counts = na.omit(counts)

####  LOOP THROUGH BATCH AND MODEL ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
batch1 <- c('ctrl1', 'ctrl2', 'dox1', 'dox2')
batch2 <- c('ctrl3', 'ctrl4', 'ctrl5', 'dox3', 'dox4', 'dox5')
for (batch in c(batch1, batch2)) {
  sample_counts1 = counts[grepl(paste(batch2, collapse='|'), names(counts))]
  for (cellline in c('SUM159', 'MCF10A', 'MDAMB231')) {
    # cellline = 'SUM159'
    sample_counts = sample_counts1[grepl(cellline, names(sample_counts1))]
    
    Rep = factor(gsub(".*_dox", "", gsub(".*_ctrl", "", names(sample_counts))))
    sample.table <- data.frame(Exp = factor(gsub('[[:digit:]]+', '', gsub(".*_", "", names(sample_counts)))),
                               Rep = factor(Rep),
                               row.names = colnames(sample_counts))
    
    ##  TRANSFORM TO DESEQ OBJECT
    sample_counts_DDS <- DESeqDataSetFromMatrix(countData =  as.matrix(round(sample_counts)),
                                         colData = sample.table, 
                                         design = ~ Exp)
    
    ##  REMOVE ZERO counts
    counts_DDS <- counts_DDS[rowSums(counts(counts_DDS)) >= 10,]
    
    ####  DIFFERENTIAL GENE EXPRESSION ####
    counts_DDS <- estimateSizeFactors(counts_DDS)
    counts_DDS.out <- DESeq(counts_DDS)
    
    counts_DDS <- estimateDispersionsGeneEst(counts_DDS)
    dispersions(counts_DDS) <- mcols(counts_DDS)$dispGeneEst
    
    counts_DDS.results <- results(counts_DDS.out)
    counts_DDS.results <- counts_DDS.results[order(counts_DDS.results$pvalue),]
    
    counts_DDS.results.volcano <- data.frame(counts_DDS.results)
    
    counts_DDS.results.volcano$GeneSymbol <- 
      human_all_genes[match(gsub("\\..*","", rownames(counts_DDS.results.volcano)),
                            human_all_genes$ensembl_gene_id), "external_gene_name"]
    
    counts_DDS.results.volcano <- counts_DDS.results.volcano[!is.na(counts_DDS.results.volcano$padj),]
    write.csv(counts_DDS.results.volcano, paste0("/Users/k2367592/Desktop/Ligand_Data/results-S3/figures/Batch2_", cellline, "_DGEX.csv"), row.names=TRUE)
  }
}

# Subset by batch
batch1 <- c('ctrl1', 'ctrl2', 'dox1', 'dox2')
batch2 <- c('ctrl3', 'ctrl4', 'ctrl5', 'dox3', 'dox4', 'dox5')
counts_1 = counts[grepl(paste(batch1, collapse='|'), names(counts))]
counts_2 = counts[grepl(paste(batch2, collapse='|'), names(counts))]


# Subset by cell line
# cellline = 'MCF10A' #### 'SUM159' # 'MDAMB231' # 'MCF10A'
counts_MCF = counts[grepl('MCF10A', names(counts))]
counts_SUM = counts[grepl('SUM159', names(counts))]
counts_MDA = counts[grepl('MDAMB231', names(counts))]



counts <- counts_1
counts <- counts_SUM



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####  PREPARE DESEQ OBJECT  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Rep = factor(gsub(".*_dox", "", gsub(".*_ctrl", "", names(counts))))
sample.table <- data.frame(Model = factor(gsub("_.*", "", names(counts))),
                           Exp = factor(gsub('[[:digit:]]+', '', gsub(".*_", "", names(counts)))),
                           Rep = factor(Rep),
                           Batch = factor(gsub("3", "2", gsub("4", "2", gsub("5", "2", gsub("2", "1", Rep))))),
                           condition = factor(gsub('.{1}$', '', names(counts))),
                           Model_batch = factor(gsub("3", "2", gsub("4", "2", gsub("5", "2", gsub("2", "1", gsub("dox", "", gsub("ctrl", "", names(counts)))))))),
                           row.names = colnames(counts))
     
##  TRANSFORM TO DESEQ OBJECT
counts_DDS <- DESeqDataSetFromMatrix(countData =  as.matrix(round(counts)),
                                     colData = sample.table, 
                                     design = ~ Exp)

##  REMOVE ZERO counts
counts_DDS <- counts_DDS[rowSums(counts(counts_DDS)) >= 10,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

####  DIFFERENTIAL GENE EXPRESSION ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##  ESTIMATE SIZE AND PERFORM DIFFERENTIAL EXPRESSION
counts_DDS <- estimateSizeFactors(counts_DDS)
counts_DDS.out <- DESeq(counts_DDS)

counts_DDS <- estimateDispersionsGeneEst(counts_DDS)
dispersions(counts_DDS) <- mcols(counts_DDS)$dispGeneEst

counts_DDS.results <- results(counts_DDS.out)
counts_DDS.results <- counts_DDS.results[order(counts_DDS.results$pvalue),]


#~~~~~~~~~~~~~~~~~~~~~~#
####  VOLCANO PLOT  ####
#~~~~~~~~~~~~~~~~~~~~~~#

counts_DDS.results.volcano <- data.frame(counts_DDS.results)

# Annotate genes
counts_DDS.results.volcano$GeneSymbol <- 
  human_all_genes[match(gsub("\\..*","", rownames(counts_DDS.results.volcano)),
                        human_all_genes$ensembl_gene_id), "external_gene_name"]

counts_DDS.results.volcano <- counts_DDS.results.volcano[!is.na(counts_DDS.results.volcano$padj),]
write.table(counts_DDS.results.volcano, file = "/Users/k2367592/Desktop/Ligand_Data/results-S3/figures/Batch1_SUM_DGEX.txt",
            sep = "\t", quote = FALSE)

counts_DDS.results.volcano$DGE <- ifelse(counts_DDS.results.volcano$log2FoldChange > 0 & counts_DDS.results.volcano$padj < 0.01, "Up",
                                         ifelse(counts_DDS.results.volcano$log2FoldChange < 0 & counts_DDS.results.volcano$padj < 0.01, "Down",
                                                "NS"))

# pdf("/Users/jelmar/Desktop/RNA_seq/results/SUM159/Volcano_SUM159.pdf")
ggplot(counts_DDS.results.volcano, aes(x = log2FoldChange, y = -log10(padj), 
                                       colour = DGE)) +
  geom_point(size = 1) +
  scale_color_manual(values=c("firebrick3", "grey", "green4")) +
  ggtitle("no BC MDA") +
  xlab("Fold change") + ylab("-log10(Q value)") +
  theme(legend.position="none")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####     TPM Counts Matrix     ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

TPM <- convertCounts(
  data.matrix(counts),
  'TPM',
  geneLength$Length
)

TPM <- as.data.frame(TPM)
TPM$GeneSymbol <- 
  human_all_genes[match(gsub("\\..*","", rownames(TPM)),human_all_genes$ensembl_gene_id), "external_gene_name"]

TPM
write.csv(TPM, "/Users/k2367592/Desktop/Ligand_Data/results-S3/figures/TPM_raw_counts.csv", row.names=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####  Remove Batch Effect  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
counts_DDS.vst <- vst(counts_DDS, blind = TRUE)

# mat <- limma::removeBatchEffect(mat, batch=counts_DDS.vst$Batch, design=mm) ### Just Batch 1 & 2
# Batch by Model (6 batches)
mat <- assay(counts_DDS.vst)
mm <- model.matrix(~condition, colData(counts_DDS.vst))
mat <- limma::removeBatchEffect(mat, batch=counts_DDS.vst$Model_batch, design=mm)
assay(counts_DDS.vst) <- mat

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####          TPM          ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

TPM <- convertCounts(
  data.matrix(assay(counts_DDS.vst)),
  'TPM',
  geneLength[rownames(mat),]
)
TPM <- as.data.frame(TPM)
TPM$GeneSymbol <- 
  human_all_genes[match(gsub("\\..*","", rownames(TPM)),human_all_genes$ensembl_gene_id), "external_gene_name"]

write.csv(TPM, "/Users/k2367592/Desktop/Ligand_Data/results-S3/figures/TPM_norm_bc_counts.csv", row.names=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####  Principal component analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
GEMM_All_PCA <- plotPCA(counts_DDS.vst, intgroup=c("Exp"), ntop = 1000, returnData = TRUE)
# pdf(paste("figures/", batch, '/', cellline, "_corrected_PCA.pdf", sep = ""))
ggplot(GEMM_All_PCA, aes(x = PC1, y = PC2, colour = group)) +
  geom_point() +
  geom_label_repel(label = row.names(GEMM_All_PCA),
                   size = 4,  force = 10, segment.size = 0.1)

# Each sample
# samples = colnames(assay(counts_DDS.vst))[grepl('SUM', colnames(assay(counts_DDS.vst)))]
# GEMM_All_PCA <- plotPCA(counts_DDS.vst[, samples], intgroup=c("Exp"), ntop = 1000, returnData = TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############# NK ligands ##############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- data.frame(
  name = factor(colnames(assay(counts_DDS.vst))),
  group = factor(gsub('_.*', '', colnames(assay(counts_DDS.vst)))),
  ctrl = factor(grepl("ctrl", colnames(assay(counts_DDS.vst)))),
  HO1 = assay(counts_DDS.vst)['ENSG00000143452.16',],
  MICA = assay(counts_DDS.vst)['ENSG00000204520.14',],
  MICB = assay(counts_DDS.vst)['ENSG00000204516.10',],
  ULBP1 = assay(counts_DDS.vst)['ENSG00000111981.5',],
  ULBP2 = assay(counts_DDS.vst)['ENSG00000131015.5',],
  ULBP3 = assay(counts_DDS.vst)['ENSG00000131019.11',],
  RAET1E = assay(counts_DDS.vst)['ENSG00000164520.11',],
  RAET1G = assay(counts_DDS.vst)['ENSG00000203722.8',],
  RAET1L = assay(counts_DDS.vst)['ENSG00000155918.8',],
  MDA5 = assay(counts_DDS.vst)['ENSG00000115267.9',],
  RIGI = assay(counts_DDS.vst)['ENSG00000107201.11',]
)

data = data[data$group == 'SUM159',]
data = data[data$group == 'MCF10A',]
data = data[data$group == 'MDAMB231',]
data$group

# row.names(assay(counts_DDS.vst))[grepl('ENSG00000107201', row.names(assay(counts_DDS.vst)))]


barplot(data$HO1, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
# abline(h=8.47, col="red") # 6.3 , 8.47
# abline(h=7.42, col="blue") # 4.55, 7.42
barplot(data$MICA, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$MICB, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$ULBP1, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$ULBP2, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$ULBP3, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$RAET1E, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$RAET1G, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$RAET1L, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$MDA5, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$RIGI, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)

#####################################

len = length(row.names((data)))
ligand <- rep( c('HO1', 'MICA', 'MICB', 'ULBP1', 'ULBP2', 'ULBP3', 'RAET1E', 'RAET1G', 'RAET1L', 'MDA5', 'RIGI'), each=len)
treatment <- rep(data$ctrl, 11)
group <- factor(data$group)
GEX <- c(data$HO1, data$MICA, data$MICB, data$ULBP1, data$ULBP2, data$ULBP3, data$RAET1E, data$RAET1G, data$RAET1L, data$MDA5, data$RIGI)
data2 = data.frame(group, ligand, treatment, GEX)
data2$ligand <- factor(data2$ligand, c('HO1', 'MICA', 'MICB', 'ULBP1', 'ULBP2', 'ULBP3', 'RAET1E', 'RAET1G', 'RAET1L', 'MDA5', 'RIGI'))
data2$treatment <- factor(data2$treatment, c('TRUE', 'FALSE'))

par(mar=c(3,4,3,1))
myplot <- boxplot(GEX ~ treatment*ligand , data=data2, 
                  boxwex=0.4 , ylab="Gene expression",
                  main="", 
                  col=c("slateblue1" , "tomato"),  
                  xaxt="n")


my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]

axis(1, 
     at = seq(1.5, 22, 2), 
     labels = my_names,
     tick=TRUE , cex.axis=1.2)

for(i in seq(0.5 , 24 , 2)){ 
  abline(v=i,lty=1, col="grey")
}
legend("bottomright", legend = c("Ctrl", "Dox"), 
       col=c("slateblue1" , "tomato"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.04, 0))






#########
genes <- c('HO1') #, 'MICA', 'MICB', 'ULBP1', 'ULBP2', 'ULBP3', 'RAET1E', 'RAET1G', 'RAET1L', 'MDA5', 'RIGI')
genes <- c('MICA', 'MICB', 'ULBP1', 'ULBP2', 'ULBP3', 'RAET1E', 'RAET1G', 'RAET1L') # p threshold = 0.002
genes <- c('MDA5', 'RIGI') # p threshold = 0.008
for (gene in genes) {
  
  par(mfrow = c(1, 1))
  expression = data[[gene]]
  cellline = data$group
  
  data$col = gsub("TRUE", "tomato", gsub("FALSE", "slateblue1", data$ctrl))
  data$pch = as.integer(gsub("TRUE", 12, gsub("FALSE", 18, data$ctrl)))
  
  plot.default(expression ~ cellline, type = "p",
               xlim = range(as.integer(unique(data$group))) + c(-0.4, 0.4), 
               xaxt = "n", col = data$col, pch = data$pch, ylab = gene)
  
  axis(1, at = seq_along(levels(data$group)), labels = levels(data$group))
  
  legend("bottomright", legend = c("Ctrl", "Dox"), 
         col=c("tomato" , "slateblue1"),
         pch = 15, bty = "n", pt.cex = 3, cex = 0.8,  horiz = F, inset = c(0, 0))
  
  ps <- c()
  
  for (group in levels(data$group)) {
    
  dat = data[c('ctrl', gene)][data$group == group,]
  colnames(dat)[colnames(dat) == gene] <- 'expression'
  dat
  res_aov <- aov(expression ~ ctrl,
                 data = dat
  )
  p = summary(res_aov)[[1]][["Pr(>F)"]][1]
  print(round(p, 5))
  if (p<0.02) {
      ps <- c(ps, "*")
    } else {
      ps <- c(ps, "ns")
    }
#    if (p<0.0005) {
#    ps <- c(ps, "***")
#  } else if (p<0.005) {
#    ps <- c(ps, "**")
#  } else if (p<0.05) {
#    ps <- c(ps, "*")
#  } else {
#    ps <- c(ps, "ns")
#  } #ps <- c(ps, round(p, 5))
  
  }
  
  axis(3, at = seq_along(levels(data$group)), labels = ps, tick = FALSE)
}


#### ALL samples, ligand comparison
head(data2)

par(mfrow = c(1, 1))
expression = data2$GEX
ligand = data2$ligand


data$col = gsub("TRUE", "tomato", gsub("FALSE", "slateblue1", data$ctrl))
data$pch = as.integer(gsub("TRUE", 12, gsub("FALSE", 18, data$ctrl)))

plot.default(expression ~ ligand, type = "p",
             xlim = range(as.integer(unique(data2$ligand))) + c(-0.4, 0.4), 
             xaxt = "n", col = data$col, pch = data$pch, ylab = 'expression')

axis(1, at = seq_along(levels(data2$ligand)), labels = levels(data2$ligand))

legend("bottomright", legend = c("Ctrl", "Dox"), 
       col=c("tomato" , "slateblue1"),
       pch = 15, bty = "n", pt.cex = 3, cex = 0.8,  horiz = F, inset = c(0, 0))

res_aov <- aov(GEX ~ treatment:group + treatment:ligand + group:ligand,
               data = data2)
summary(res_aov)

TukeyHSD(res_aov, which = "treatment:group")
TukeyHSD(res_aov, which = "treatment:ligand")
TukeyHSD(res_aov, which = "group:ligand")

res_aov <- aov(GEX ~ ligand:treatment,
               data = data2)
summary(res_aov)
TukeyHSD(res_aov, which = "ligand:treatment")



### REGRESSION ###
library("ggpubr")

ligand = 'RAET1L'

plot.mat = matrix(c(1, 1,
                    2, 3),
                  nrow = 2, byrow = T)
layout(plot.mat)
par(mar = c(2,2,2,2))
cor.test(data$HO1, data[[ligand]], method=c("pearson", "kendall", "spearman"))
plot(data$HO1, data[[ligand]])
barplot(data$HO1, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)
barplot(data$MICA, col=data$ctrl, legend = TRUE, names.arg=data$name, las=2)







dev.off()




res.aov <- aov(ULBP3 ~ ctrl, data = data)
summary(res.aov)
TukeyHSD(res.aov)


assay(counts_DDS.vst)['ENSG00000143452.16',]












#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####  CHECK HORMAD1 expression ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

HOcounts = assay(counts_DDS.vst)['ENSG00000143452.16',]
plot(density(HOcounts, bw = 0.4))

library(diptest)
library(LaplacesDemon)
dip.test(HOcounts)
is.unimodal(HOcounts)
is.bimodal(HOcounts)
is.trimodal(HOcounts)



