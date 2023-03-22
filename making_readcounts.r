setwd ("/genetica_2/excerpt_files/4_batchs/align_samples")

#chamar tabela com nome de todas as amostras ex 10928_W1
table <- read.table("ids.txt", header = FALSE) 

#pegar todos os readcounts de miRNAs maduros
txt_files_ls <- list.files( ".", pattern="*miRNAmature_sense.txt", all.files=FALSE, full.names=FALSE, recursive = TRUE)
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T)})

txtfilenames <- table$V1

#juntar as colunas diferentes (IDs de participantes) com as mesmas linhas (iD de miRNA)
for (i in seq_along(txt_files_df)){
  colnames(txt_files_df[[i]]) <- c("ReferenceID", paste(txtfilenames[i],"uni",sep = ""), paste(txtfilenames[i],"TotalCount",sep = ""),
                                   paste(txtfilenames[i],"multi",sep = ""), paste(txtfilenames[i],"nada",sep = ""))
}

#tabelas com read counts
library(plyr)
a <- join_all(txt_files_df, by = "ReferenceID",type = "full")

#tabela total counts
b <- a[,grepl(pattern = "*TotalCount", x = names(a))|grepl(pattern = "ReferenceID", x = names(a))]
b[is.na(b)] <- 0
write.table(b, "4b_ReadCountTotais.txt", col.names = T, row.names = F, quote = F)

#tabela com uniquely mapped counts
c  <- a[,grepl(pattern = "*uni", x = names(a))|grepl(pattern = "ReferenceID", x = names(a))]
c[is.na(c)] <- 0
write.table(c, "4b_ReadCountUniquely.txt", col.names = T, row.names = F, quote = F)

#tabela com adjusted multimapped counts 
d  <- a[,grepl(pattern = "*multi", x = names(a))|grepl(pattern = "ReferenceID", x = names(a))]
d[is.na(d)] <- 0
write.table(d, "4b_ReadCountMultiadjusted.txt", col.names = T, row.names = F, quote = F)



######################
setwd("/genetica_1/Jessica/mestrado/resultados_seq")
unicount <- read.csv("4b_ReadCountUniquely.txt", row.names=1, sep="")
#######################
# CHAMAR BIBLIOTECAS #
#######################
library("DESeq2")
library ("ggplot2")
library("ggfortify")
library('EnhancedVolcano')
library("ggrepel")
library("dplyr")



#################################
# PREPARACAO DAS TABELAS INPUT #
#################################
# Retirada das colunas do batch 1
uni_todos <- select (unicount, -8)
uni_inpd <- select (unicount, -8, -25, -26, -27, -28)

#Preparacao da tabela de readcounts (countData)
mRC_uni <- as.matrix(uni_inpd)
mRC_uni <- as.data.frame (mRC_uni)

colnames(mRC_uni) <- gsub("uni", "", colnames(mRC_uni))
colnames(mRC_uni) <- gsub("X", "", colnames(mRC_uni))

sampleinfo <- read.csv2("sample_information_140521.csv", h=T)
#arrumar o que é factor, numeric etc
#deixar o id igual ao do dSEQ
sampleinfo$ID <- gsub("_", ".", sampleinfo$ID)

#deixar a sampleinfo na mesma ordem que o seq
colnames(mRC_uni) -> ids
as.data.frame(ids) -> ids
colnames(ids) <- "ID"
join(ids, sampleinfo, by="ID", type="inner") -> sampleinfo_order
rownames(sampleinfo_order) <- sampleinfo_order$ID
#pra nao confundir
rm(sampleinfo)
#checar se a ordem está igual
rownames(sampleinfo_order)
colnames(mRC_uni)

#tirar a primeira coluna de id
sampleinfo_todosgroup <- sampleinfo_order[, c(-2,-3)]


####################################
# ANLISE  DE EXPRESSAO DIFERENCIAL #
####################################
#DESeq2 quick start -> design da analise = todos analise grupo
dds_uni <- DESeqDataSetFromMatrix(countData = mRC_uni,
                                  colData = sampleinfo_todosgroup,
                                  design = ~ batch + sex + age + group)

#Selecionar apenas miRNAs com pelo menos 24 amostras com no minimo 3 reads 
keep_counts <- rowSums(counts(dds_uni) >= 3) >= 18
dds_uni <- dds_uni[keep_counts,]

#Extrair a matriz de reads normalizada
dds_uni <- estimateSizeFactors(dds_uni)
norm_mRC_uni <- counts(dds_uni,normalized = T)

#Renomear as linhas e colunas pra ficar mais bonito
rownames(norm_mRC_uni) <- gsub("[:MIMAT].*", "", rownames(norm_mRC_uni))

ntd_uni <- normTransform(dds_uni)


######## PCS #######
pdf("pca_normtransform.pdf")
plotPCA(ntd_uni, intgroup="group")
plotPCA(ntd_uni, intgroup="sex")
plotPCA(ntd_uni, intgroup="batch")
plotPCA(ntd_uni, intgroup="wave")
dev.off()

#parece que o mais apropriado é o rlog
pdf("plotpca_rlog.pdf")
plotPCA(rld, intgroup="group")
plotPCA(rld, intgroup="sex")
plotPCA(rld, intgroup="batch")
plotPCA(rld, intgroup="wave")
dev.off()
#By default the function uses the top 500 most variable genes. You can change this by adding the ntop argument and specifying how many genes you want to use to draw the plot.

#NOTE: The plotPCA() function will only return the values for PC1 and PC2. If you would like to explore the additional PCs in your data or if you would like to identify genes that contribute most to the PCs, you can use the prcomp() function. For example, to plot any of the PCs we could run the following code:
  
# Input is a matrix of log transformed values
rld <- rlog(dds_uni, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
#grupo
df <- cbind(sampleinfo_order, pca$x)
PC1toPC5_var<- pca$sdev*pca$sdev
PC1toPC5_per <- round(PC1toPC5_var/sum(PC1toPC5_var)*100,1)
pc1 <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, col = group), size = 3) +xlab(paste("PC1 - ", PC1toPC5_per[1], "%", sep = "")) + ylab(paste("PC2 - ", PC1toPC5_per[2], "%", sep = ""))
pc2 <- ggplot(df) + geom_point(aes(x=PC2, y=PC3, col = group), size = 3)+ xlab(paste("PC2 - ", PC1toPC5_per[2], "%", sep = "")) + ylab(paste("PC3 - ", PC1toPC5_per[3], "%", sep = ""))
library(ggpubr)
pc3 <- ggplot(df) + geom_point(aes(x=PC3, y=PC4, col = group), size = 3) +xlab(paste("PC3 - ", PC1toPC5_per[3], "%", sep = "")) + ylab(paste("PC4 - ", PC1toPC5_per[4], "%", sep = ""))
pc4 <- ggplot(df) + geom_point(aes(x=PC4, y=PC5, col = group), size = 3) +xlab(paste("PC4 - ", PC1toPC5_per[4], "%", sep = "")) + ylab(paste("PC5 - ", PC1toPC5_per[5], "%", sep = ""))

pdf(file="pca_rlog.pdf")
ggarrange(pc1, pc2, pc3, pc4,common.legend = TRUE,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
dev.off()

##Heatmap##
library("pheatmap")
select_uni <- order(rowMeans(counts(dds_uni,normalized=TRUE)),
                    decreasing=TRUE)[1:50]

df_uni <- as.data.frame(colData(dds_uni)[,c("Sex","condition")])

rownames(df_uni) <- gsub("[:MIMAT].*", "", rownames(df_uni))
rownames(ntd_uni) <- gsub("[:MIMAT].*", "", rownames(ntd_uni))
colnames(ntd_uni) <- gsub("[:Total].*","", colnames(ntd_uni))

pheatmap(assay(ntd_uni)[select_uni,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df_uni, fontsize_row=8, fontsize_col
         =8, main="Heatmap - top 50 expressed miRNA")


#####################
### GRUPOS NA W1 ###
#####################

#Definir "Control" como ref para comparacao
dds_uni$condition <- relevel(dds_uni$condition, ref = "A")
#results(dds, c("condition", "A", "C"))

# lists the coefficients
dds_uni <- DESeq(dds_uni)
resultsNames(dds_uni)
BA <- results(dds_uni, name="group_B_vs_A")
CA <- results(dds_uni, name="group_C_vs_A")
DA <- results(dds_uni, name="group_D_vs_A")
sex <- results(dds_uni, name="sex_2_vs_1")
batch <- results(dds_uni, name="batch_2_vs_1")
age <- results(dds_uni, name="age")

age
summary(age)
BA
CA
DA

#############################
### COMPARACAO D x A #####
#############################
group_uni <- as.data.frame(DA)
summary(group_uni)

##VOLCANO PLOT##
#Criar coluna com nome de miRNAs
group_uni$miRNA <- rownames(group_uni)
group_uni$miRNA <- sub("[:MIMAT].*", "",group_uni$miRNA)

#Retirada de NAs
sgroup_uni <- group_uni[complete.cases(group_uni), ]

#Cria��o de coluna de up ou down regulated
sgroup_uni$diffexpressed <- "No difference"
sgroup_uni$diffexpressed[sgroup_uni$log2FoldChange > 0.6 & sgroup_uni$padj < 0.05] <- "Up regulated"
sgroup_uni$diffexpressed[sgroup_uni$log2FoldChange < -0.6 & sgroup_uni$padj < 0.05] <- "Down regulated"

sgroup_uni$delabel <- NA
sgroup_uni$delabel[sgroup_uni$diffexpressed != "No difference"] <- sgroup_uni$miRNA[sgroup_uni$diffexpressed != "No difference"]

#Gr�fico
pdf("volcanogroup.pdf")
ggplot(sgroup_uni, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point(aes(color = diffexpressed)) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-0.6, 0.6), col="gray") +
  geom_hline(yintercept=-log10(0.05), col="gray") + xlim (-5,5) + 
  scale_color_manual(values= c("cyan4", "black", "coral1")) +
  ggtitle("Differential Expressed miRNAs: Group comparison") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Dif. Expressed miRNA') + 
  geom_text_repel()
dev.off()

#####################
### GRUPOS NA W2 ###
#####################


#################################
### COMPARACAO W1 x W2 #####
#################################
#Ajuste da tabela de read counts (OPCIONAL tirar possiveis outliers e seus pares
mRCt_uni <- mRC_uni [,-7]
colnames(mRCt_uni) -> ids
as.data.frame(ids) -> ids
colnames(ids) <- "ID"
sampleinfo <- read.csv2("sample_information_140521.csv", h=T)
join(ids, sampleinfo, by="ID", type="inner") -> sampleinfo_order
rownames(sampleinfo_order) <- sampleinfo_order$ID
#pra nao confundir
rm(sampleinfo)
#checar se a ordem está igual
rownames(sampleinfo_order)
colnames(mRCt_uni)

#tirar colunas que não pertencem nessa analise. Deixar: wave, batch, sex, age, subjectid (entender se deixa o grupo tb)
sampleinfo_paired <- sampleinfo_order[, c(-1,-2)]

#DESeq2 quick start -> design da analise
ddst_uni <- DESeqDataSetFromMatrix(countData = mRCt_uni,
                                   colData = sampleinfo_paired,
                                   design = ~ subjectid + wave)

#Selecionar apenas miRNAs com pelo menos 10 amostras com no minimo 10 reads 
keept_uni <- rowSums(counts(ddst_uni) >= 3) >= 18
ddst_uni <- ddst_uni[keept_uni,]

#Extrair a matriz de reads normalizada
ddst_uni <- estimateSizeFactors(ddst_uni)
norm_mRCt_uni <- counts(ddst_uni,normalized = T)

#Renomear as linhas pra ficar mais bonito
rownames(norm_mRCt_uni) <- gsub("[:MIMAT].*", "", rownames(norm_mRCt_uni))

#Definir "Control" como ref para comparacao
ddst_uni$wave <- relevel(ddst_uni$wave, ref = "w1")

# lists the coefficients
ddst_uni <- DESeq(ddst_uni)
resultsNames(ddst_uni)

### Analise diferencial #####
wave_uni <- results(ddst_uni, name="wave_w2_vs_w1")
wave_uni <- as.data.frame(wave_uni)
summary(wave_uni)

##VOLCANO PLOT##
#Criar coluna com nome de miRNAs
wave_uni$miRNA <- rownames(wave_uni)
wave_uni$miRNA <- sub("[:MIMAT].*", "",wave_uni$miRNA)

#Retirada de NAs
swave_uni <- wave_uni[complete.cases(wave_uni), ]

#Cria��o de coluna de up ou down regulated
swave_uni$diffexpressed <- "No difference"
swave_uni$diffexpressed[swave_uni$log2FoldChange > 0.6 & swave_uni$pvalue < 0.05] <- "Up regulated"
swave_uni$diffexpressed[swave_uni$log2FoldChange < -0.6 & swave_uni$pvalue < 0.05] <- "Down regulated"

swave_uni$delabel <- NA
swave_uni$delabel[swave_uni$diffexpressed != "No difference"] <- swave_uni$miRNA[swave_uni$diffexpressed != "No difference"]

#Gr�fico --> falta adicionar label aos demiRNAs
pdf("wave_volcano_notadj.pdf")
ggplot(swave_uni, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(aes(color = diffexpressed)) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-0,6, 0.6), col="gray") +
  geom_hline(yintercept=-log10(0.05), col="gray") + xlim (-10,10) + 
  scale_color_manual(values= c("cyan4","black","coral1" )) +
  ggtitle("Differential Expressed miRNAs: w2 x w1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Dif. Expressed miRNA') + 
  geom_text_repel()
dev.off()



## SALVAR OS RESULTADOS EM CSV ##
write.csv(HC_FEP_uni, "HC_FEP_uni.csv", col.names=TRUE, row.names = TRUE)
write.csv(HC_FEP8W_uni, "HC_FEP8W_uni.csv", col.names=TRUE, row.names = TRUE)
write.csv(wave_uni, "wave_uni.csv", col.names=TRUE, row.names = TRUE)

############################################################
############################################################
############################################################

## GRAFICOS
#PCA -> usando o valor dos reads absolutos
rownames(sdesign_uni) <- gsub("[:Total].*","", rownames(sdesign_uni))

pca_condition_uni <- counts(dds_uni, normalized=T)
pca_cond_uni <- prcomp(t(pca_condition_uni))

## para calcular a % que o PC representa:
pca_cond.var_uni <- pca_cond_uni$sdev^2
pca_cond.var.per_uni <- round(pca_cond.var_uni/sum(pca_cond.var_uni)*100,1)

### PCA - Batch effect ###
pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,2])
pdf("pca1x2_batch_uninormalized.pdf")
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sampleinfo_todosgroup$batch)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC1 x PC2, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,2], Y=pca_cond_uni$x[,3])
pdf("pca2x3_batch_uninormalized.pdf")
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sampleinfo_todosgroup$batch)) +
  geom_text() +
  xlab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC2 x PC3, coloured by ext Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,3], Y=pca_cond_uni$x[,4])
pdf("pca3x4_batch_uninormalized.pdf")
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sampleinfo_todosgroup$batch)) +
  geom_text() +
  xlab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_cond.var.per_uni[4], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC3 x PC4, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,3])
pdf("pca1x3_batch_uninormalized.pdf")
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sampleinfo_todosgroup$batch)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC1 x PC3, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

### PCA - SEXO ###
pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,2])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign$Sex)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC1 x PC2, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,2], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign_uni$Sex)) +
  geom_text() +
  xlab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC2 x PC3, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,3], Y=pca_cond_uni$x[,4])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign_uni$Sex)) +
  geom_text() +
  xlab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_cond.var.per_uni[4], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC3 x PC4, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign_uni$Sex)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC1 x PC3, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

### PCA - CONDITION ###
pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,2])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign$condition)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Condition") +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC1 x PC2, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,2], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign_uni$condition)) +
  geom_text() +
  xlab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC2 x PC3, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,3], Y=pca_cond_uni$x[,4])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign_uni$condition)) +
  geom_text() +
  xlab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_cond.var.per_uni[4], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Condition") +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC3 x PC4, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = sdesign_uni$condition)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Condition") +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC1 x PC3, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

#####################################
##Violin Plot - miRNAs Dif. expr.##
#miRNA 30d-5p#
norm_mRCt_uni <- as.data.frame (t(norm_mRC_uni))
colnames(norm_mRCt_uni) <- gsub("[-]", "_", colnames(norm_mRCt_uni))
Group = c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W" )
mir <- norm_mRCt_uni$hsa_miR_30d_5p
vp <- cbind( "miRNA" = mir, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)

ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA") + xlab ("Grupo") + ggtitle("miRNA mir-30d-5p") +
  theme(plot.title = element_text(hjust = 0.5))

#miRNA let7b-5p#
norm_mRCt_uni <- as.data.frame (t(norm_mRC_uni))
colnames(norm_mRCt_uni) <- gsub("[-]", "_", colnames(norm_mRCt_uni))
Group = c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W","FEP-8W" )
mir <- norm_mRCt_uni$hsa_let_7b_5p
vp <- cbind( "miRNA" = mir, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)

ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA") + xlab ("Grupo") + ggtitle("miRNA let 7b-5p") +
  theme(plot.title = element_text(hjust = 0.5))

#### VAN###
library(ggpubr)
# para selecionar o gene mais significante:
topGene <- rownames(results)[which.min(results$padj)]
plot <- plotCounts(dds, gene=topGene, intgroup=c("group","wave", "CBCL_level", "subjectID"), 
                   returnData=TRUE)

# Retirar o grupo NotClassified (n�o interessa para o grafico)
plot %>% filter(group != "NotClassified") %>% ggline(x = "wave", y = "count", add = c("violin", "mean_se"),color = "group", facet.by = "group", palette = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF")) +
  ylab("Normalized counts") +
  ggtitle(topGene, "NR3C1, Incident vs Control, pvalue = 0.004, padj = 0.345, log2FC=-0.352") + theme_gray()


#####################################
#####################################

## PREDICAO DE ALVOS DOS miRNAs - mirnatap ##

#carregar bibliotecas
library(miRNAtap)
library(topGO)
library(org.Hs.eg.db)

#predicted target analysis
#Predicao de alvos#
mir = 'miR-30d-5p'
predictions = getPredictedTargets(mir, species = 'hsa', method = 'geom', min_src = 2)
head(predictions)

#An�lise de vias#
rankedGenes = predictions[,'rank_product']
selection = function(x) TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez")
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes,annot = annFUN.GO2genes, GO2genes = allGO2genes,geneSel = selection, nodeSize=10)
results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks

allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 40)
allRes[,c('GO.ID','Term','KS')]
write.csv(allRes, "path_anal_uni.csv", row.names = TRUE)

#Transformar EntrezID em Gene Symbol#
hs <- org.Hs.eg.db
my.symbols <- row.names(predictions)
geneid <- select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
row.names(predictions) <- geneid$SYMBOL
lista_genes <- row.names(predictions)
write.csv(lista_genes, "predicted_targets.txt")