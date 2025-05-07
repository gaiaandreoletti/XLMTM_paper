## This script generates the CPM data matrix, performs unsupervised (PCA and unsupervised clustering) and identifies 
## differentially expressed genes in all comparisons between samples of cohort #2.
## Unsupervised analysis is performed retaining the highly variable genes, i.e. those gene with a coefficient of variation
## larger than the "percent" percentile of the coefficient of variation distribution.
## Supervised analysis if performed using edgeR glmQLFTest with a paired design on all possible comparisons between time points.
## Differentially expressed genes are defined as genes with an absolute fold-change larger than FC.thr and a q-value smaller
## than qvalue.thr.
## This script also creates the data, annotation, and class files for GSEA.
## Results (files and plots) are saved into the directory "../results/cohort2".

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
  qvalue.thr<-0.05
  FC.thr<-2

### directories and libraries ###
  sub.column<-"cohort2"
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

### sub-sample info file and data matrix ###
  infos<-sample.infos[sample.infos[,sub.column]=="yes",]
  count.data<-data.count.matrix[,infos$Samplename]
  setwd(res.dir)

  ### filter and normalize ###
    group.column<-"time"
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
    autoplot(pca,x=1,y=2,data=infos,colour='time',scale=0,
             size=4,label=TRUE,label.size=3,label.hjust=-0.2)
    ggsave(paste0("Unsupervised_",percent,"_PCA.pdf"),useDingbats=FALSE,width=8,height=7,units="in",dpi=600)
    ### Clustering and heatmap ###
    annot.col<-data.frame(Time=infos$time,Tissue=infos$biopsy_tissue,Mutation=infos$mutation_type)
    rownames(annot.col)<-colnames(filter.data)
    #show_col(hue_pal()(40))
    Time<-hue_pal()(3)
    names(Time)<-levels(annot.col$Time)
    Tissue<-brewer.pal(n=11,name='Paired')[c(3,9,11)]
    names(Tissue)<-levels(annot.col$Tissue)
    Mutation<-brewer.pal(n=9,name='YlGnBu')[c(1,9)]
    names(Mutation)<-levels(annot.col$Mutation)
    ann.colors<-list(Time=Time,Tissue=Tissue,Mutation=Mutation)
    #heat.colors<-colorRampPalette(c("blue", "white", "red"))(256)
    pheatmap(filter.data,scale="row",clustering_distance_rows="correlation",
             clustering_distance_cols="correlation",clustering_method="average",
             cluster_rows=T, cluster_cols=T,border_color="grey80",treeheight_row=0,
             annotation_col=annot.col,annotation_colors=ann.colors,annotation_names_col=F,
             labels_col=infos$Samplename,labels_row=rep("",dim(filter.data)[1]),
             filename = paste0("Unsupervised_",percent,"_clustering.pdf"),width=8,height=7)
    dev.off() 
  
  ### Compute DEGs ###
    setwd(sup.dir)
    comparisons<-combinations(3,2,unique(infos[,group.column]))
    for (i in 1:dim(comparisons)[1]){
      time1<-comparisons[i,1]
      time2<-comparisons[i,2]
      sub.infos<-infos[infos$time==time1|infos$time==time2,]
      patient<-factor(sub.infos$patient)
      treatment.time<-factor(sub.infos$time)
      design<-model.matrix(~patient+treatment.time)
      sub.DGE.obj<-DGE.obj[,sub.infos$Samplename]
      sub.DGE.obj<-estimateDisp(sub.DGE.obj,design,robust=TRUE)
      sqrt(sub.DGE.obj$common.dispersion)
      plotBCV(sub.DGE.obj)
      fit<-glmQLFit(sub.DGE.obj,design) 
      DEGs.obj<-glmQLFTest(fit)
      degs<-topTags(DEGs.obj,n=dim(DEGs.obj$table)[1])
      summary(decideTests(DEGs.obj))
      plotMD(DEGs.obj)
      abline(h=c(-1,1),col="blue")
      degs$table$FC<-2^degs$table$logFC
      degs$table$FC<-2^degs$table$logFC 
      degs$table$FC<-ifelse(degs$table$FC>1,degs$table$FC,-1/degs$table$FC)
      degs<-data.frame(rownames(DGE.obj),degs$table[rownames(DGE.obj),c(7,2,5:6)])
      colnames(degs)[1]<-"GeneSymbol"
      length(rownames(degs[degs$FC>=2&degs$FDR<=0.05,]))
      length(rownames(degs[degs$FC<=-2&degs$FDR<=0.05,]))
      write.table(degs,file=paste0(sub.column,"_",time2,"-vs-",time1,".txt",sep=""),sep="\t", row.names = F)
      # Clustering on DEGs ##
      deg.genes<-rownames(degs[(degs$FC<=-FC.thr|degs$FC>=FC.thr)&degs$FDR<=qvalue.thr,])
      ## clustering compared samples ##
      data.to.clus<-cpm.log.matrix[deg.genes,sub.infos$Samplename]
      sub.annot.col<-data.frame(Time=sub.infos$time,Tissue=sub.infos$biopsy_tissue,Mutation=sub.infos$mutation_type)
      rownames(sub.annot.col)<-colnames(data.to.clus)
      sub.ann.colors<-list(Time=Time[unique(sub.infos$time)],
                       Tissue=Tissue[unique(sub.infos$biopsy_tissue)],
                       Mutation=Mutation[unique(sub.infos$mutation_type)])
      if (dim(data.to.clus)[1]<20){font.size=10}
      if (dim(data.to.clus)[1]>20&dim(data.to.clus)[1]<90){font.size=5}
      if (dim(data.to.clus)[1]>90){font.size=1}
      heat1<-pheatmap(data.to.clus,scale="row",clustering_distance_rows="correlation",
               clustering_distance_cols="correlation",clustering_method="average",
               cluster_rows=T, cluster_cols=T,border_color="grey80",treeheight_row=0,
               annotation_col=sub.annot.col,annotation_colors=sub.ann.colors,annotation_names_col=F,
               labels_col=sub.infos$Samplename,cutree_rows =2,cutree_cols=2,fontsize_row=font.size)
       ## clustering with all samples ##
       s1<-grep(unique(infos[,group.column])[1],colnames(cpm.log.matrix))
       s2<-grep(unique(infos[,group.column])[2],colnames(cpm.log.matrix))
       s3<-grep(unique(infos[,group.column])[3],colnames(cpm.log.matrix))
       gap1<-length(s1)
       gap2<-gap1+length(s2)
       gap3<-gap2+length(s3)
       data.to.clus<-cpm.log.matrix[deg.genes,c(s1,s2,s3)]
       if (dim(data.to.clus)[1]<20){
          heat2<-pheatmap(data.to.clus,scale="row",clustering_distance_rows="correlation",
                clustering_distance_cols="correlation",clustering_method="average",
                cluster_rows=T, cluster_cols=F,border_color="grey80",treeheight_row=0,
                annotation_col=annot.col,annotation_colors=ann.colors,annotation_names_col=F,
                labels_col=colnames(data.to.clus),cutree_rows=2,gaps_col=c(gap1,gap2,gap3),
                cellwidth=15,cellheight=15,fontsize_row=font.size)
       }else{
          heat2<-pheatmap(data.to.clus,scale="row",clustering_distance_rows="correlation",
                 clustering_distance_cols="correlation",clustering_method="average",
                 cluster_rows=T, cluster_cols=F,border_color="grey80",treeheight_row=0,
                 annotation_col=annot.col,annotation_colors=ann.colors,annotation_names_col=F,
                 labels_col=colnames(data.to.clus),cutree_rows=2,gaps_col=c(gap1,gap2,gap3),fontsize_row=font.size)
       }
       plot_list=list()
       plot_list[[1]] = heat1[[4]]
       plot_list[[2]] = heat2[[4]]
       g<-grid.arrange(arrangeGrob(grobs=plot_list,ncol=2))
       ggsave(file=paste0(sub.column,"_DEG_",time2,"-vs-",time1,".pdf"),g,width=12,height=7,useDingbats=FALSE)
      dev.off()
    }

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
    
