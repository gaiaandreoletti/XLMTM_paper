setwd("/fast/analysis/20210518_asipiro_normal_skeletal/analysis_ventilation/")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("tximport")
library("tximport")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/slow/data/moved_from_fast/20210302_Nantes_RNAseq/gencode.v37.annotation.gtf")
k <- keys(txdb, keytype = "GENEID")
df <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
head(tx2gene)

a<- gsub("\\..*","",tx2gene[,1])
b<- gsub("\\..*","",tx2gene[,2])
c<-cbind(a,b)
colnames(c)=colnames(tx2gene)
tx2gene <- as.data.frame(c)
head(tx2gene)

dir= getwd()


#BiocManager::install("DESeq2")
library(DESeq2)

# data_counts_cases <- read.table(file = "star_deseq_raw_counts_batch6.txt", header = TRUE, row.names = 1)
 data_counts_cases <- read.table(file = "../ASPIRO_raw_count.txt", header = TRUE, row.names = 1)
dim(data_counts_cases)
head(data_counts_cases)
# data_counts_cases$Gene_a <- gsub("\\..*","",rownames(data_counts_cases))
# library(biobroom)
# library(dplyr)
# library(annotables)
# data_counts_cases__annotation <- data_counts_cases %>% inner_join(grch38, by=c("Gene_a"="ensgene")) 
# dim(data_counts_cases__annotation)
# rownames(data_counts_cases__annotation) <- make.unique(data_counts_cases__annotation$symbol)
# head(data_counts_cases__annotation)
# data_counts_cases_symbol <- data_counts_cases__annotation[,1:45]
# head(data_counts_cases_symbol)
# 
# colnames(data_counts_cases_symbol) <- gsub("F([0-9, A-Z]+)_.*", "\\1", colnames(data_counts_cases_symbol))
# head(data_counts_cases_symbol)

 data_counts_cases_baseline <- data_counts_cases[,c("F3300074001",
                                                 "F330008J001",
                                                "F33000B4001",
                                                 "F330001Y001",
                                                  "F33000J5001", #out
                                               "F33000N3001",
                                                 "F33000W0001",
                                                "F33001JV001",
                                                 "F33000LT001",
                                                 "F330021N001",
                                                 "F330027P001",
                                                 "F33001RX001", #out??
                                                 "F33001RY001",
                                                 "F330027Q001", #out
                                                 "F33004QG001")]


# data_counts_cases_baseline <- data_counts_cases_symbol[,c("F33000LT001",
#                                                           "F33000B4001",
#                                                           "F330008J001",
#                                                           "F3300074001",
#                                                           "F330027P001",
#                                                           "F330027Q001",
#                                                           "F33004QG001",
#                                                           "F330021N001",
#                                                           "F33000N3001",
#                                                           "F330001Y001",
#                                                           "F33000W0001",
#                                                           "F33001JV001",
#                                                           "F33001RX001",
#                                                           "F33001RY001",
#                                                           "F33000J5001",
#                                                           "F330021K001",
#                                                           "F33002ZW001")]

dim(data_counts_cases_baseline)

data_counts_cases_wk24 <- data_counts_cases[,c("F330008H001",
                                               "F33000RM001",
                                               "F33000RH001",
                                               "F33000W1001",
                                               "F33000YU001",
                                               "F33001E8001",
                                               "F33001E6001",
                                               "F33001RZ001",
                                               "F330021L001",
                                               "F33001ZV001",
                                               "F330040W001",
                                               "F33002ZZ001",
                                               "F33002ZY001",
                                               "F33003RF001",
                                               "F33004R6001")] 

dim(data_counts_cases_wk24)

data_counts_cases_wk48 <- data_counts_cases[,c("F33000RK001", 
                                               "F3300075001",
                                               "F33001DB001",
                                               "F33001E7001",
                                               "F33000YV001",
                                               "F330021M001",
                                               "F33001ZW001",
                                               "F330023D001",
                                               "F3300331001",
                                               "F33002JE001",
                                               "F33004R4001",
                                               "F3300458001",
                                               "F3300459001",
                                               "F33004R7001")] 

dim(data_counts_cases_wk48)

#baseline
samples <- read.table("ASPIRO_infos_vent.txt", header = T)
rownames(samples) <- samples$sampleID

library(data.table)
samples_baseline <- samples[samples$time %like% "baseline", ]
samples_baseline <- samples_baseline[- grep("decreased", samples_baseline$treatment),]
samples_wk48 <- samples[samples$time %like% "wk48", ]
samples_wk24 <- samples[samples$time %like% "wk24", ]

write.table(samples_baseline, "samples_baseline.txt", quote = FALSE)
write.table(samples_wk48, "samples_wk48.txt", quote = FALSE)
write.table(samples_wk24, "samples_wk24.txt", quote = FALSE)


ncol(data_counts_cases_baseline) == nrow(samples_baseline)
all_cohort <- data_counts_cases_baseline[,order(colnames(data_counts_cases_baseline))]
samples <- samples_baseline[order(rownames(samples_baseline)),]
colnames(all_cohort) == rownames(samples)

write.csv(all_cohort,"all_cohort_counts_baseline.csv", quote = FALSE)

#wk24
ncol(data_counts_cases_wk24) == nrow(samples_wk24)
all_cohort_wk24 <- data_counts_cases_wk24[,order(colnames(data_counts_cases_wk24))]
samples_wk24 <- samples_wk24[order(rownames(samples_wk24)),]
colnames(all_cohort_wk24) == rownames(samples_wk24)

write.csv(all_cohort_wk24,"all_cohort_counts_all_wk24.csv", quote = FALSE)

#wk48
ncol(data_counts_cases_wk48) == nrow(samples_wk48)
samples_wk48 <- samples_wk48[- grep("decreased", samples_wk48$treatment),]
all_cohort_wk48 <- data_counts_cases_wk48[,order(colnames(data_counts_cases_wk48))]
samples_wk48 <- samples_wk48[order(rownames(samples_wk48)),]
colnames(all_cohort_wk48) == rownames(samples_wk48)

write.csv(all_cohort_wk48,"all_cohort_counts_all_wk48.csv", quote = FALSE)


######### baseline ########
# samples$Replicate <- as.factor(samples$Replicate)
samples_baseline$Responsiveness_simplified <- as.factor(samples_baseline$Responsiveness_simplified)

##############
#calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(all_cohort)

#pull PCA values out of the PCA object
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

#assign labels to the data frame
pca_uncorrected[,"condition"] = samples$Condition
pca_uncorrected[,"library_method"] = samples$batch
pca_uncorrected[,"sampleID"] = samples$sampleID

#plot the PCA
#create a classic 2-dimension PCA plot (first two principal components) with conditions and library methods indicated
cols <- c("Controls" = "#481567FF", "Aspiro_baseline" = "#1F968BFF")
p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA, RNA-seq counts uncorrected data", color="Condition", shape="Library Method")
p1 = p1 + scale_colour_manual(values = cols)

## batch correction combact 
# BiocManager::install("BatchQC")
## combact ###
# devtools::install_github("zhangyuqing/sva-devel")
library(sva)
adjusted <- ComBat_seq(all_cohort, batch=samples$batch, group=NULL)
dim(adjusted)
# corrected_data = ComBat_seq(counts = as.matrix(all_cohort), batch = samples$Condition, group = groups)
# corrected_data = cbind(uncorrected_data[,c("Gene","Chr")], corrected_data)
#calculate principal components for the uncorrected data
pca_corrected_obj = prcomp(adjusted[,samples$sampleID])

#pull PCA values out of the PCA object
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)

#assign labels to the data frame
pca_corrected[,"Condition"] = samples$Condition
pca_corrected[,"batch"] = samples$batch
pca_corrected[,"name"] = samples$sampleID

#as above, create a PCA plot for comparison to the uncorrected data
cols <- c("Controls" = "#481567FF", "Aspiro_baseline" = "#1F968BFF")
p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=Condition, shape=batch, label=name))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq counts for Controls/Aspiro_baseline samples (batch corrected data Combact)", color="Condition", shape="batch")
 p2 = p2 + scale_colour_manual(values = cols) #+ geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
p2
pdf(file="Uncorrected-vs-BatchCorrected-PCA.pdf")
library(gridExtra)
grid.arrange(p1, p2, nrow = 2)
dev.off()

write.csv(all_cohort,"all_cohort_counts_all_adjusted.csv", quote = FALSE)
##############

## remove outlier in baseline F330027Q001 & F330021N001

all_cohort <- select(all_cohort, -contains(c("F330027Q001","F330021N001")))
samples_baseline <- samples_baseline[!grepl("F330027Q001", samples_baseline$sampleID),]
samples_baseline <- samples_baseline[!grepl("F330021N001", samples_baseline$sampleID),]

dds <- DESeqDataSetFromMatrix(countData = round(all_cohort),
                              colData = samples_baseline,
                              design = ~  Responsiveness_simplified)


cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
summary(dds)

dds <- estimateSizeFactors(dds)
nsamp <- 12 #90% of 15 
keep <- rowSums(counts(dds, normalized=TRUE)> 5)>= nsamp 
table(keep) # 16870 
dds <- dds[keep,]  
dds <- DESeq(dds, fitType='local')
summary(dds)
resultsNames(dds) 
#res1 <- results(dds)
res1 <- results(dds, contrast = c("Responsiveness_simplified","unresponsive","responsive"))
head(res1,2)
dim(res1)
table(res1$padj<0.05)#42   #8
res1 <- res1[order(res1$padj), ] 
head(res1,2)
res1_data <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res1_data)[1] <- "Gene" 
head(res1_data,2)
table(res1$padj<0.05 & res1$log2FoldChange > 1)# 26  up in unresponsive - down in responsive #3
table(res1$padj<0.05 & res1$log2FoldChange < - 1)# 3 up in responsive  #3
library(biobroom)
library(dplyr)
library(annotables)
res_tidy <- tidy.DESeqResults(res1)
res_tidy__annotation <- res_tidy %>% inner_join(grch38, by=c("gene"="ensgene")) 
date_is <-Sys.time()
write.csv(res_tidy__annotation, file=paste(date_is,"diffexpr-results_Responsiveness_simplified_unresponsive_vs_responsive_OutliersRemoved_baseline.csv"))
vsd <- vst(dds, blind=FALSE, fitType='local')

library(ggrepel)

pdf("PCA2_baseline.pdf", paper = "USr")
# DESeq2::plotPCA(vsd, intgroup="Replicate")+ ggtitle("Replicate") 
# DESeq2::plotPCA(vsd, intgroup="sample_ID")+ ggtitle("sample_ID")
DESeq2::plotPCA(vsd, intgroup="Responsive")+ ggtitle("Responsive")
DESeq2::plotPCA(vsd, intgroup="hours_on_ventilation")+ ggtitle("hours_on_ventilation")
DESeq2::plotPCA(vsd, intgroup="treatment")+ ggtitle("treatment") +   geom_text_repel(aes(label = name), max.overlaps = Inf) +
  geom_point(color = 'red') +
  theme_classic(base_size = 16)
dev.off()
# 
# assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$Condition)
# DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition post correction") + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
# DESeq2::plotPCA(vsd, intgroup="sampleID")+ ggtitle("sampleID")
# dev.off()


# rlog <- assay(rlog(dds,blind=T))
# mod <- model.matrix(~ Condition, samples)
# cbat <- ComBat(rlog, batch=cond$Condition, mod=mod)

library(pheatmap)
pdf("heatmap_baseline.pdf")
padj.cutoff <- 0.05 
lfc.cutoff <- 2
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res1$threshold <- threshold     
sigOE <- data.frame(subset(res1, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res1, threshold==TRUE))
df <- as.data.frame(colData(dds)[,c("Responsive")])
rownames(df)<- colData(dds)$sampleID
pheatmap(assay(vsd)[select,], show_rownames=T, annotation_col=df, cutree_cols = 2, cutree_rows = 2, fontsize = 12, fontsize_row  = 12, scale = "row",width = 2, height = 1)
dev.off()

# d<-assay(vsd)[select,names]
# write.table(d,"table.txt", quote = FALSE)
# d <- read.table("table.txt", header = T, row.names = 1)
# pheatmap(d, show_rownames=T,  fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf("Volcano_baseline.pdf")	  
EnhancedVolcano(res1, lab = rownames(res1), 
                x = 'log2FoldChange',      y = 'padj', 
                xlim = c(-13, 13),      pCutoff = 0.05, 
                FCcutoff = 2) 
dev.off()

#pdf("Enrichment_CD4_racewhite.pdf",width=20,height=8)
library(clusterProfiler)
library(org.Hs.eg.db)
OrgDb <- org.Hs.eg.db # can also be other organisms 
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res1$symbol = mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="ENSEMBL", 
                     keytype="SYMBOL", 
                     multiVals="first")
res1$entrez = mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="ENTREZID", 
                     keytype="SYMBOL", 
                     multiVals="first")
res1$name =   mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="GENENAME", 
                     keytype="SYMBOL", 
                     multiVals="first")

# DESeq results to pathways in 60 Seconds with the fgsea package ###########
library(fgsea)
ens2symbol = mapIds(org.Hs.eg.db, 
                    keys=res1_data$Gene,  
                    column="ENSEMBL", 
                    keytype="SYMBOL", 
                    multiVals="first")

ens2symbol <- tibble::enframe(ens2symbol)
ens2symbol
res <- inner_join(res1_data, ens2symbol, by=c("Gene"="name"))

ranks <- res$log2FoldChange
names(ranks) <- res$Gene
head(ranks, 20)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("/slow/data/moved_from_fast/20210302_Nantes_RNAseq/h.all.v7.0.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up/down reg pathways unresponsive vs responsive baseline") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

fwrite(topPathways, file ="myDT_baseline.csv")

############ week 24 #########
# samples_wk24$Replicate <- as.factor(samples_wk24$Replicate)

ncol(data_counts_cases_wk24) == nrow(samples_wk24)
all_cohort_wk24 <- data_counts_cases_wk24[,order(colnames(data_counts_cases_wk24))]
samples_wk24 <- samples_wk24[order(rownames(samples_wk24)),]
colnames(all_cohort_wk24) == rownames(samples_wk24)

# all_cohort_wk24 <- select(all_cohort_wk24, -contains("F33001ZV001"))
# samples_wk24 <- samples_wk24[!grepl("F33001ZV001", samples_wk24$sampleID),]


#all_cohort_wk24
dds <- DESeqDataSetFromMatrix(countData = round(all_cohort_wk24),
                              colData = samples_wk24,
                              design = ~  Responsiveness_simplified)

cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
summary(dds)

dds <- estimateSizeFactors(dds)
nsamp <- 14 #90% of 15 
keep <- rowSums(counts(dds, normalized=TRUE)> 3)>= nsamp 
table(keep) # 16757 
dds <- dds[keep,]  
dds <- DESeq(dds, fitType='local')
summary(dds)
resultsNames(dds) 
#res1 <- results(dds)
res1 <- results(dds, contrast = c("Responsiveness_simplified","unresponsive","responsive"))
head(res1,2)
dim(res1)
table(res1$padj<0.05)# 1070 
res1 <- res1[order(res1$padj), ] 
head(res1,2)
res1_data <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res1_data)[1] <- "Gene" 
head(res1_data,2)
table(res1$padj<0.05 & res1$log2FoldChange > 1)# 687   up in unresponsive - down in responsive #3
table(res1$padj<0.05 & res1$log2FoldChange < - 1)# 241  up in responsive  
library(biobroom)
library(dplyr)
library(annotables)
res_tidy <- tidy.DESeqResults(res1)
res_tidy__annotation <- res_tidy %>% inner_join(grch38, by=c("gene"="ensgene")) 
date_is <-Sys.time()
write.csv(res_tidy__annotation, file=paste(date_is,"diffexpr-results_Responsiveness_simplified_unresponsive_vs_responsive_wk24.csv"))
vsd <- vst(dds, blind=FALSE, fitType='local')

library(ggrepel)

pdf("PCA2_wk24.pdf", paper = "USr")
# DESeq2::plotPCA(vsd, intgroup="Replicate")+ ggtitle("Replicate") 
# DESeq2::plotPCA(vsd, intgroup="sample_ID")+ ggtitle("sample_ID")
DESeq2::plotPCA(vsd, intgroup="Responsive")+ ggtitle("Responsive")
DESeq2::plotPCA(vsd, intgroup="hours_on_ventilation")+ ggtitle("hours_on_ventilation")
DESeq2::plotPCA(vsd, intgroup="treatment")+ ggtitle("treatment") +   geom_text_repel(aes(label = name), max.overlaps = Inf) +
  geom_point(color = 'red') +
  theme_classic(base_size = 16)
dev.off()
# 
# assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$Condition)
# DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition post correction") + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
# DESeq2::plotPCA(vsd, intgroup="sampleID")+ ggtitle("sampleID")
# dev.off()


# rlog <- assay(rlog(dds,blind=T))
# mod <- model.matrix(~ Condition, samples)
# cbat <- ComBat(rlog, batch=cond$Condition, mod=mod)

library(pheatmap)
pdf("heatmap_wk24.pdf")
padj.cutoff <- 0.05 
lfc.cutoff <- 2
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res1$threshold <- threshold     
sigOE <- data.frame(subset(res1, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res1, threshold==TRUE))
df <- as.data.frame(colData(dds)[,c("Responsive")])
rownames(df)<- colData(dds)$sampleID
pheatmap(assay(vsd)[select,], show_rownames=T, annotation_col=df, cutree_cols = 2, cutree_rows = 3, fontsize = 4, fontsize_row  = 1, scale = "row",width = 2, height = 1)
dev.off()


a <- read.table("week24_subset.txt")
subset_genes <- a$V1
pheatmap(assay(vsd)[subset_genes,], show_rownames=T, annotation_col=df, cutree_cols = 2, cutree_rows = 2, fontsize = 4, fontsize_row  = 5, scale = "row",width = 2, height = 1)



# d<-assay(vsd)[select,names]
# write.table(d,"table.txt", quote = FALSE)
# d <- read.table("table.txt", header = T, row.names = 1)
# pheatmap(d, show_rownames=T,  fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf("Volcano_wk24.pdf")	  
EnhancedVolcano(res1, lab = rownames(res1), 
                x = 'log2FoldChange',      y = 'padj', 
                xlim = c(-13, 13),      pCutoff = 0.05, 
                FCcutoff = 2) 
dev.off()

#pdf("Enrichment_CD4_racewhite.pdf",width=20,height=8)
library(clusterProfiler)
library(org.Hs.eg.db)
OrgDb <- org.Hs.eg.db # can also be other organisms 
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res1$symbol = mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="ENSEMBL", 
                     keytype="SYMBOL", 
                     multiVals="first")
res1$entrez = mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="ENTREZID", 
                     keytype="SYMBOL", 
                     multiVals="first")
res1$name =   mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="GENENAME", 
                     keytype="SYMBOL", 
                     multiVals="first")

# DESeq results to pathways in 60 Seconds with the fgsea package ###########
library(fgsea)
ens2symbol = mapIds(org.Hs.eg.db, 
                    keys=res1_data$Gene,  
                    column="ENSEMBL", 
                    keytype="SYMBOL", 
                    multiVals="first")

ens2symbol <- tibble::enframe(ens2symbol)
ens2symbol
res <- inner_join(res1_data, ens2symbol, by=c("Gene"="name"))

ranks <- res$log2FoldChange
names(ranks) <- res$Gene
head(ranks, 20)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("/slow/data/moved_from_fast/20210302_Nantes_RNAseq/h.all.v7.0.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up/down reg pathways unresponsive vs responsive wk24") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

fwrite(topPathways, file ="myDT_wk24.csv")

####
########## week 48 #########

ncol(data_counts_cases_wk48) == nrow(samples_wk48)
all_cohort_wk48 <- data_counts_cases_wk48[,order(colnames(data_counts_cases_wk48))]
samples_wk48 <- samples_wk48[order(rownames(samples_wk48)),]
colnames(all_cohort_wk48) == rownames(samples_wk48)


 all_cohort_wk48 <- select(all_cohort_wk48, -contains("F330023D001"))
 samples_wk48 <- samples_wk48[!grepl("F330023D001", samples_wk48$sampleID),]

dds <- DESeqDataSetFromMatrix(countData = round(all_cohort_wk48),
                              colData = samples_wk48,
                              design = ~  Responsiveness_simplified)


cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
summary(dds)

dds <- estimateSizeFactors(dds)
nsamp <- 13 #90% of 15 
keep <- rowSums(counts(dds, normalized=TRUE)> 5)>= nsamp 
table(keep) # 16870 
dds <- dds[keep,]  
dds <- DESeq(dds, fitType='local')
summary(dds)
resultsNames(dds) 
#res1 <- results(dds)
res1 <- results(dds, contrast = c("Responsiveness_simplified","unresponsive","responsive"))
head(res1,2)
dim(res1)
table(res1$padj<0.05)#84   
res1 <- res1[order(res1$padj), ] 
head(res1,2)
res1_data <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res1_data)[1] <- "Gene" 
head(res1_data,2)
table(res1$padj<0.05 & res1$log2FoldChange > 1)# 38  up in unresponsive - down in responsive 
table(res1$padj<0.05 & res1$log2FoldChange < - 1)# 44 up in responsive  
library(biobroom)
library(dplyr)
library(annotables)
res_tidy <- tidy.DESeqResults(res1)
res_tidy__annotation <- res_tidy %>% inner_join(grch38, by=c("gene"="ensgene")) 
date_is <-Sys.time()
write.csv(res_tidy__annotation, file=paste(date_is,"diffexpr-results_Responsiveness_simplified_unresponsive_vs_responsive_OutliersRemoved_wk48.csv"))
vsd <- vst(dds, blind=FALSE, fitType='local')

library(ggrepel)

pdf("PCA_wk48.pdf", paper = "USr")
# DESeq2::plotPCA(vsd, intgroup="Replicate")+ ggtitle("Replicate") 
# DESeq2::plotPCA(vsd, intgroup="sample_ID")+ ggtitle("sample_ID")
DESeq2::plotPCA(vsd, intgroup="Responsive")+ ggtitle("Responsive")
DESeq2::plotPCA(vsd, intgroup="hours_on_ventilation")+ ggtitle("hours_on_ventilation")
DESeq2::plotPCA(vsd, intgroup="treatment")+ ggtitle("treatment") +   geom_text_repel(aes(label = name), max.overlaps = Inf) +
  geom_point(color = 'red') +
  theme_classic(base_size = 16)
dev.off()
# 
# assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$Condition)
# DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition post correction") + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
# DESeq2::plotPCA(vsd, intgroup="sampleID")+ ggtitle("sampleID")
# dev.off()


# rlog <- assay(rlog(dds,blind=T))
# mod <- model.matrix(~ Condition, samples)
# cbat <- ComBat(rlog, batch=cond$Condition, mod=mod)

library(pheatmap)
pdf("heatmap_wk48.pdf")
padj.cutoff <- 0.05 
lfc.cutoff <- 2
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res1$threshold <- threshold     
sigOE <- data.frame(subset(res1, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res1, threshold==TRUE))
df <- as.data.frame(colData(dds)[,c("Responsive")])
rownames(df)<- colData(dds)$sampleID
pheatmap(assay(vsd)[select,], show_rownames=T, annotation_col=df, cutree_cols = 2, cutree_rows = 2, fontsize = 12, fontsize_row  = 8, scale = "row",width = 2, height = 1)
dev.off()

# d<-assay(vsd)[select,names]
# write.table(d,"table.txt", quote = FALSE)
# d <- read.table("table.txt", header = T, row.names = 1)
# pheatmap(d, show_rownames=T,  fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf("Volcano_wk48.pdf")	  
EnhancedVolcano(res1, lab = rownames(res1), 
                x = 'log2FoldChange',      y = 'padj', 
                xlim = c(-13, 13),      pCutoff = 0.05, 
                FCcutoff = 2) 
dev.off()

# DESeq results to pathways in 60 Seconds with the fgsea package ###########
library(fgsea)
ens2symbol = mapIds(org.Hs.eg.db, 
                    keys=res1_data$Gene,  
                    column="ENSEMBL", 
                    keytype="SYMBOL", 
                    multiVals="first")

ens2symbol <- tibble::enframe(ens2symbol)
ens2symbol
res <- inner_join(res1_data, ens2symbol, by=c("Gene"="name"))

ranks <- res$log2FoldChange
names(ranks) <- res$Gene
head(ranks, 20)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("/slow/data/moved_from_fast/20210302_Nantes_RNAseq/h.all.v7.0.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up/down reg pathways unresponsive vs responsive wk48") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

fwrite(topPathways, file ="myDT_wk48.csv")

###
dds <- estimateSizeFactors(dds)
nsamp <- 22 #90% of 25
keep <- rowSums(counts(dds, normalized=TRUE)> 10)>= nsamp 
table(keep) # 13002
dds <- dds[keep,]  
dds <- DESeq(dds, fitType='local')
summary(dds)
resultsNames(dds) 
#res1 <- results(dds)
res1 <- results(dds, contrast = c("Condition","Controls","wk48"))
head(res1,2)
dim(res1)
table(res1$padj<0.05)#7054   
res1 <- res1[order(res1$padj), ] 
head(res1,2)
res1_data <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res1_data)[1] <- "Gene" 
head(res1_data,2)
table(res1$padj<0.05 & res1$log2FoldChange > 2)# 446   up in WT - down in Aspiro
table(res1$padj<0.05 & res1$log2FoldChange < - 2)# 324 up in Aspiro 
library(biobroom)
library(dplyr)
library(annotables)
res_tidy <- tidy.DESeqResults(res1)
res_tidy__annotation <- res_tidy %>% inner_join(grch38, by=c("gene"="ensgene")) 
date_is <-Sys.time()
write.csv(res_tidy__annotation, file=paste(date_is,"diffexpr-results_Controls_vs_Aspiro_week48.csv"))
vsd <- vst(dds, blind=FALSE, fitType='local')
pdf("PCA2_week48.pdf")
# DESeq2::plotPCA(vsd, intgroup="Replicate")+ ggtitle("Replicate") 
# DESeq2::plotPCA(vsd, intgroup="sample_ID")+ ggtitle("sample_ID")
DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition pre correction")
assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$Condition)
DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition post correction") + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
DESeq2::plotPCA(vsd, intgroup="sampleID")+ ggtitle("sampleID")
dev.off()


# rlog <- assay(rlog(dds,blind=T))
# mod <- model.matrix(~ Condition, samples)
# cbat <- ComBat(rlog, batch=cond$Condition, mod=mod)

library(pheatmap)
pdf("heatmap_week48.pdf")
padj.cutoff <- 0.05 
lfc.cutoff <- 2
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res1$threshold <- threshold     
sigOE <- data.frame(subset(res1, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res1, threshold==TRUE))
df <- as.data.frame(colData(dds)[,c("Condition")])
rownames(df)<- colData(dds)$sampleID
pheatmap(assay(vsd)[select,], show_rownames=F, annotation_col=df, cutree_cols = 2, cutree_rows = 2, fontsize = 12, fontsize_row  = 12, scale = "row",width = 2, height = 1)
dev.off()

# d<-assay(vsd)[select,names]
# write.table(d,"table.txt", quote = FALSE)
# d <- read.table("table.txt", header = T, row.names = 1)
# pheatmap(d, show_rownames=T,  fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf("Volcano_week48.pdf")	  
EnhancedVolcano(res1, lab = rownames(res1), 
                x = 'log2FoldChange',      y = 'padj', 
                xlim = c(-13, 13),      pCutoff = 0.05, 
                FCcutoff = 2) 
dev.off()

#pdf("Enrichment_CD4_racewhite.pdf",width=20,height=8)
library(clusterProfiler)
library(org.Hs.eg.db)
OrgDb <- org.Hs.eg.db # can also be other organisms 
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res1$symbol = mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="ENSEMBL", 
                     keytype="SYMBOL", 
                     multiVals="first")
res1$entrez = mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="ENTREZID", 
                     keytype="SYMBOL", 
                     multiVals="first")
res1$name =   mapIds(org.Hs.eg.db, 
                     keys=row.names(res1),  
                     column="GENENAME", 
                     keytype="SYMBOL", 
                     multiVals="first")

# DESeq results to pathways in 60 Seconds with the fgsea package ###########
library(fgsea)
ens2symbol = mapIds(org.Hs.eg.db, 
                    keys=res1_data$Gene,  
                    column="ENSEMBL", 
                    keytype="SYMBOL", 
                    multiVals="first")

ens2symbol <- tibble::enframe(ens2symbol)
ens2symbol
res <- inner_join(res1_data, ens2symbol, by=c("Gene"="name"))

ranks <- res$log2FoldChange
names(ranks) <- res$Gene
head(ranks, 20)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("/fast/analysis/20210302_Nantes_RNAseq/h.all.v7.0.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up/down reg pathways Controls vs Aspiro week48") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

fwrite(topPathways, file ="myDT_week48.csv")

plotGseaTable(pathways.hallmark[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)



# rownames(samples) <- colnames(txi.salmon$counts)
# samples$Conidtion <- as.factor(samples$Conidtion)
# colnames(all_480_counts) == rownames(colData)
# counts <- txi.salmon$counts
# dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~Conidtion)

### Differential Splicing Analysis ####
BiocManager::install("DEXSeq")
install.packages("dplyr")
library("DEXSeq")
BiocManager::install("rnaseqDTU")
library("rnaseqDTU")
samps <- samples
head(samps)
samps$Conidtion <- factor(samps$Conidtion)
table(samps$Conidtion)
files <- file.path(dir, "quants", samples$run_accession_salmon, "quant.sf")
# files <- file.path("/path/to/dir", samps$sample_id, "quant.sf")
names(files) <- samps$sample_ID
head(files)
library(tximport)
txi <- tximport(files, type="salmon", txOut=TRUE,
                countsFromAbundance="scaledTPM")
# cts <- txi$counts
cts <- txi.salmon$counts
cts <- cts[rowSums(cts) > 0,]
cts <- cts[,-c(1:6)]

# library(GenomicFeatures)
# gtf <- "gencode.v37.annotation.gtf"
# txdb.filename <- "gencode.v37.annotation.sqlite"
# txdb <- makeTxDbFromGFF(gtf)
# saveDb(txdb, txdb.filename)
# 
# txdb <- loadDb(txdb.filename)
# txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
# tab <- table(txdf$GENEID)
# txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# write.csv(tx2gene,"tx2gene.csv", quote = FALSE)

cts[1:3,1:3]
range(colSums(cts)/1e6) # 27.64240 49.36784 million paired-end reads that were mapped to the transcriptome using Salmon.
all(rownames(cts) %in% tx2gene$GENEID)
txdf <- tx2gene[match(rownames(cts),tx2gene$GENEID),]
# txdf <- na.exclude(txdf)
all(rownames(cts) == txdf$GENEID)

counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)

library(DRIMSeq)
samps <- samps[-c(1:6),]
d <- dmDSdata(counts=counts, samples=samps)
d
methods(class=class(d))
counts(d[1,])[,1:4]

n <- 0
n.small <- 0
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=0,
              min_samps_feature_prop=n.small, min_feature_prop=0,
              min_samps_gene_expr=n, min_gene_expr=0)
d
table(table(counts(d)$gene_id))
design_full <- model.matrix(~Conidtion, data=DRIMSeq::samples(d))
colnames(design_full)

set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="ConidtionMT_U7")
})

res <- DRIMSeq::results(d)
head(res)

res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)

###
library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))

write.table(cts,"countFiles.txt",quote = FALSE, sep = "\t")
countFiles = list.files(pattern=".quant$", full.names=TRUE)
countFiles <- files

#python /home/gandreoletti/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_prepare_annotation.py Homo_sapiens.GRCh38.103.gtf Homo_sapiens.GRCh38.DEXSeq.gff
flattenedFile = list.files(pattern="gtf$", full.names=TRUE)
flattenedFile <- flattenedFile[2]

dxd = DEXSeqDataSet(
  round(cts),
  sampleData=samps,
  design= ~ exon + Conidtion,
  featureID = txdf$GENEID,
  groupID = txdf$TXNAME)

split( seq_len(ncol(dxd)), colData(dxd)$exon )
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dxd = estimateSizeFactors( dxd, geoMeans=geoMeans )
dxd = estimateDispersions( dxd,fitType="local" )
mplotDispEsts( dxd )

system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~Condition + exon)
})

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE, normalized=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)
