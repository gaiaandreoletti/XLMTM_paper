## This script generates the CPM data matrix, performs unsupervised (PCA and unsupervised clustering) and identifies 
## differentially expressed genes in the comparison between baseline samples of cohort #1 and of cohort #2.
## Unsupervised analysis is performed retaining the highly variable genes, i.e. those gene with a coefficient of variation
## larger than the "percent" percentile of the coefficient of variation distribution.
## Supervised analysis if performed using edgeR exactTest.
## Differentially expressed genes are defined as genes with an absolute fold-change larger than FC.thr and a q-value smaller
## than qvalue.thr.
## This script also creates the data, annotation, and class files for GSEA.
## Results (files and plots) are saved into the directory "../results/baselines".

## local.dir: change based on the directory where data are stored
## project: if the project name is changed, then the file name of annotation, sample info, and data have to be changed accordingly
## percent: the percentile of the CV distribution to define the highly variable genes in unsupervised analysis
## FC.thr: the fold change absolute threshold to define the differentially expressed genes
## qvalue.thr: the q-value threshold to define the differentially expressed genes

### R-3.6.2
### edgeR_3.28.1

#### Settings ####
## general path settings
  local.dir<-"your local directory"
  project<-"ASPIRO"
  main.path<-file.path(local.dir,project)
## settings for the unsupervised analysis
  percent<-0.95
## settings for plotting differential expressed genes
  FC.thr<-2
  qvalue.thr<-0.05

### directories and libraries ###
  sub.column<-"baselines"
  data.dir<-file.path(main.path,"raw_data")
  main.res.dir<-file.path(main.path,"results")
  if (!file.exists(main.res.dir)){dir.create(main.res.dir)}
  res.dir<-file.path(main.res.dir,sub.column)
  gsea.dir<-file.path(res.dir,"GSEA")
  unsup.dir<-file.path(res.dir,"unsupervised")
  sup.dir<-file.path(res.dir,"supervised")
  if (!file.exists(res.dir)){dir.create(res.dir)}
  if (!file.exists(gsea.dir)){dir.create(gsea.dir)}
  if (!file.exists(unsup.dir)){dir.create(unsup.dir)}
  if (!file.exists(sup.dir)){dir.create(sup.dir)}
  library(edgeR)
  library(RColorBrewer)
  library(gplots)
  library(gtools)
  library(ggfortify)
  library(pheatmap)
  library(scales)
  library(gridExtra)

###	import raw count data ###
  raw.count.file<-paste0(project,"_raw_count.txt")
  info.file<-paste0(project,"_infos.txt")
  annot.file<-paste0(project,"_annot.txt")
  setwd(data.dir)
  data.count.matrix<-read.table(file=raw.count.file,header =TRUE,sep ="\t",row.names=1)
  sample.infos<-read.table(file=info.file,header =TRUE,sep ="\t",stringsAsFactors = F)
  annot<-read.table(file=annot.file,header =TRUE,sep ="\t")
  data.count.matrix<-data.count.matrix[,sample.infos$sampleID]
  colnames(data.count.matrix)<-sample.infos$Samplename

### info file and data matrix sub-sampling ###
  infos<-sample.infos[sample.infos[,sub.column]=="yes",]
  count.data<-data.count.matrix[,infos$Samplename]
  setwd(res.dir)

  ### filter and normalize ###
    group.column<-"dose_level"
    DGE.obj<-DGEList(counts=count.data,group=infos[,group.column],genes=rownames(count.data))
    dim(DGE.obj)
  # filter genes with less than 1 cpm in at least as many samples as in the smaller type
    min.n.samples<-min(table(infos[,group.column])) # minimal number of replicates
    DGE.obj<-DGE.obj[rowSums(cpm(DGE.obj)>1)>=min.n.samples,]  
    dim(DGE.obj)
  # re-compute the library sizes on the filtered DGElist
    DGE.obj$samples$lib.size<-colSums(DGE.obj$counts)
  # calculate normalization factors to scale the raw library sizes
    DGE.obj<-calcNormFactors(DGE.obj)
    DGE.obj$samples
  
  ### Compute counts per million (CPM) and moderated log2(CPM) ### 
    cpm.matrix<-cpm(DGE.obj)
    cpm.log.matrix<-cpm(DGE.obj,log=TRUE,prior.count=1)
    cpm.matrix<-data.frame("GeneSymbols"=rownames(cpm.matrix),cpm.matrix)
    write.table(cpm.matrix,file=paste0(sub.column,"_cpm.txt"),sep="\t",quote=F,row.names=F)
    cpm.log.matrix.to.write<-data.frame("GeneSymbols"=rownames(cpm.log.matrix),cpm.log.matrix)
    write.table(cpm.log.matrix.to.write,file=paste0(sub.column,"_log2cpm.txt"),col.names=NA,sep="\t",quote=F)
  
  ### Compute fragment per kilobase (FPKM) if pair-end or read per kilobase (RPKM) if single-end ### 
    genelength<-annot$Length
    names(genelength)<-annot$GeneID
    sum(names(genelength)%in%rownames(DGE.obj))
    genelength.reordered<-genelength[rownames(DGE.obj)]
    rpkm.matrix<-rpkm(DGE.obj,gene.length=genelength.reordered)
    rpkm.matrix<-data.frame("GeneSymbols"=rownames(rpkm.matrix),rpkm.matrix,stringsAsFactors=F)
    write.table(rpkm.matrix,file=paste0(sub.column,"_rpkm.txt"),sep="\t",quote=F,row.names=F)
  
  ### unsupervised analysis: clustering and PCA ###
    setwd(unsup.dir)
    GEdata.matrix<-cpm.log.matrix  
    cvs<-apply(GEdata.matrix,1,sd)/apply(GEdata.matrix,1,mean)
    filter.data<-GEdata.matrix[cvs>=quantile(cvs,percent),]
    ### PCA ###
    pca<-prcomp(t(filter.data),center=T,scale.=T)
    summary(pca)
    rownames(infos)<-infos$Samplename
    autoplot(pca,x=1,y=2,data=infos,colour='dose_level',scale=0,shape="dose_level",
             size=3,label=TRUE,label.size=3,label.hjust=-0.2)
    ggsave(paste0("Unsupervised_",percent,"_PCA.pdf"),useDingbats=FALSE,width=8,height=7,units="in",dpi=600)
    dev.off() 
    ### Clustering and heatmap ###
    annot.col<-data.frame(Dose=infos$dose_level,Mutation=infos$mutation_type,Age=infos$AgeLess2)
    rownames(annot.col)<-colnames(filter.data)
    #show_col(hue_pal()(40))
    Dose<-brewer.pal(n=9,name='YlOrRd')[c(4,6)]
    names(Dose)<-levels(annot.col$Dose)
    Mutation<-brewer.pal(n=9,name='YlGnBu')[c(1,9)]
    names(Mutation)<-levels(annot.col$Mutation)
    Age<-brewer.pal(n=11,name='Paired')[c(3,9)]
    names(Age)<-levels(annot.col$Age)
    ann.colors<-list(Dose=Dose,Mutation=Mutation,Age=Age)
    pheatmap(filter.data,scale="row",clustering_distance_rows="correlation",
             clustering_distance_cols="correlation",clustering_method="average",
             cluster_rows=T, cluster_cols=T,border_color="grey80",treeheight_row=0,
             annotation_col=annot.col,annotation_colors=ann.colors,annotation_names_col=F,
             labels_col=infos$Samplename,labels_row=rep("",dim(filter.data)[1]),
             filename = paste0("Unsupervised_",percent,"_clustering.pdf"),width=8,height=7)
    
  ### Compute DEGs ###
    setwd(sup.dir)
    ctrl<-"cohort1"
    treat<-"cohort2"
    DGE.obj<-estimateDisp(DGE.obj)
    plotBCV(DGE.obj)
    sqrt(DGE.obj$common.dispersion)
    DEGs.obj<-exactTest(DGE.obj, pair=c(ctrl,treat))
    degs<-topTags(DEGs.obj,n=dim(DEGs.obj$table)[1])
    degs$table$FC<-2^degs$table$logFC 
    degs$table$FC<-ifelse(degs$table$FC>1,degs$table$FC,-1/degs$table$FC)
    degs<-data.frame(rownames(DGE.obj),degs$table[rownames(DGE.obj),c(6,2,4:5)])
    colnames(degs)[1]<-"GeneSymbol"
    write.table(degs,file=paste0(sub.column,"_",treat,"-vs-",ctrl,".txt"),sep="\t", row.names = F)
  
  ## Clustering on DEGs ##
    deg.genes<-rownames(degs[(degs$FC<=-FC.thr|degs$FC>=FC.thr)&degs$FDR<=qvalue.thr,])
    data.to.clus<-cpm.log.matrix[deg.genes,infos$Samplename]
    annot.col<-data.frame(Dose=infos$dose_level,Mutation=infos$mutation_type,Age=infos$AgeLess2)
    rownames(annot.col)<-colnames(data.to.clus)
    ann.colors<-list(Dose=Dose,Mutation=Mutation,Age=Age)
    if (dim(data.to.clus)[1]<20){font.size=10}
    if (dim(data.to.clus)[1]>20&dim(data.to.clus)[1]<90){font.size=5}
    if (dim(data.to.clus)[1]>90){font.size=1}
    pheatmap(data.to.clus,scale="row",clustering_distance_rows="correlation",
                    clustering_distance_cols="correlation",clustering_method="average",
                    cluster_rows=T, cluster_cols=T,border_color="grey80",treeheight_row=0,
                    annotation_col=annot.col,annotation_colors=ann.colors,annotation_names_col=F,
                    labels_col=infos$Samplename,cutree_rows =2,cutree_cols=2,fontsize_row=font.size,
                    filename = paste0(sub.column,"_",treat,"-vs-",ctrl,".pdf"),width=8,height=7)
 
### creation of files for GSEA ###
    setwd(gsea.dir)
    # creation of the data file #
    write.table(cpm.matrix,file="data.txt",sep="\t",quote=F,row.names=F)
    # creation of .chip file #
    chipannot<-data.frame(annot$GeneID,annot$GeneID,annot$GeneID)
    colnames(chipannot)<-c("Probe Set ID","Gene Symbol","Gene Title")
    write.table(chipannot,"annot.chip",sep="\t",row.names=F,quote=F)
    # creation of .cls file #
    level.temp<-levels(as.factor(infos[,group.column]))
    cls.vct<-as.vector(infos[,group.column])
    cls.label<-level.temp[match(level.temp,unique(cls.vct))]
    filename<-paste(group.column,".cls",sep="")
    cat(c(length(cls.vct),length(cls.label),1),"\n",sep="\t",file=filename)
    cat(c("#",cls.label),"\n",sep="\t",file=filename,append=T)
    cat(cls.vct,sep="\t",file=filename,append=T)
   
