#' @title prog_DEG_selection
#' @description Prognostic DEGs selection based on survival data for subsequent process
#' @details Input takes several data frames which are the exports of MCKNNV and MCKNNV_tSNE
#' @param RNA_seq_dataframe Data frame of the RNA sequencing profiles from public database
#' @param clin_info Data frame of the clinical information containing full survival status and survival time.
#' @param MCKV_group the manually grouping results based on MCKNNV_tSNE
#' @export


prog_DEG_selection <- function(RNA_seq_dataframe, clin_info, MCKV_group){
  Type = data.frame(MCKV_group[,1])
  row.names(Type) <- row.names(MCKV_group)
  colnames(Type)[1] <- "Subtype"
  RNA_seq_dataframe=as.matrix(RNA_seq_dataframe)

  #Comparing group1 and group3
  #extracting sample information of different groups
  survLow=Type[Type$Subtype=="1",,drop=F]
  survHigh=Type[Type$Subtype=="3",,drop=F]
  RNA_seq_dataframe_Low=RNA_seq_dataframe[,row.names(survLow)]
  RNA_seq_dataframe_High=RNA_seq_dataframe[,row.names(survHigh)]
  RNA_seq_dataframe_Low_High=cbind(RNA_seq_dataframe_Low,RNA_seq_dataframe_High)
  RNA_seq_dataframe_Low_High <- as.matrix(RNA_seq_dataframe_Low_High)
  dimnames=list(rownames(RNA_seq_dataframe_Low_High),colnames(RNA_seq_dataframe_Low_High))
  RNA_seq_dataframe_Low_High=matrix(as.numeric(as.matrix(RNA_seq_dataframe_Low_High)),nrow=nrow(RNA_seq_dataframe_Low_High),dimnames=dimnames)
  RNA_seq_dataframe_Low_High=RNA_seq_dataframe_Low_High[rowMeans(RNA_seq_dataframe_Low_High)>1,]#过滤掉表达量低的基因

  conNum=ncol(RNA_seq_dataframe_Low)
  treatNum=ncol(RNA_seq_dataframe_High)
  Type_by_order=c(rep(1,conNum), rep(2,treatNum))

  #Differentially expressed genes analysis
  logFCfilter=0.585
  fdrFilter=0.05
  DEG_outTab_1vs3=data.frame()
  for(i in row.names(RNA_seq_dataframe_Low_High)){
    unigene_exp=data.frame(expression=RNA_seq_dataframe_Low_High[i,], Type=Type_by_order)
    wilcoxTest=wilcox.test(expression ~ Type, data=unigene_exp)
    pvalue=wilcoxTest$p.value
    conGeneMeans=mean(RNA_seq_dataframe_Low_High[i,1:conNum])
    treatGeneMeans=mean(RNA_seq_dataframe_Low_High[i,(conNum+1):ncol(RNA_seq_dataframe_Low_High)])
    logFC=treatGeneMeans-conGeneMeans
    conMed=median(RNA_seq_dataframe_Low_High[i,1:conNum])
    treatMed=median(RNA_seq_dataframe_Low_High[i,(conNum+1):ncol(RNA_seq_dataframe_Low_High)])
    diffMed=treatMed-conMed
    if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
      DEG_outTab_1vs3=rbind(DEG_outTab_1vs3,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
    }
  }
  pValue=DEG_outTab_1vs3[,"pValue"]
  fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
  DEG_outTab_1vs3=cbind(DEG_outTab_1vs3, fdr=fdr)
  Sig_DEG_outTab_1vs3=DEG_outTab_1vs3[( abs(as.numeric(as.vector(DEG_outTab_1vs3$logFC)))>logFCfilter & as.numeric(as.vector(DEG_outTab_1vs3$fdr))<fdrFilter),]

  #Comparing group2 and group3
  #extracting sample information of different groups
  survLow=Type[Type$Subtype=="2",,drop=F]
  survHigh=Type[Type$Subtype=="3",,drop=F]
  RNA_seq_dataframe_Low=RNA_seq_dataframe[,row.names(survLow)]
  RNA_seq_dataframe_High=RNA_seq_dataframe[,row.names(survHigh)]
  RNA_seq_dataframe_Low_High=cbind(RNA_seq_dataframe_Low,RNA_seq_dataframe_High)
  RNA_seq_dataframe_Low_High <- as.matrix(RNA_seq_dataframe_Low_High)
  dimnames=list(rownames(RNA_seq_dataframe_Low_High),colnames(RNA_seq_dataframe_Low_High))
  RNA_seq_dataframe_Low_High=matrix(as.numeric(as.matrix(RNA_seq_dataframe_Low_High)),nrow=nrow(RNA_seq_dataframe_Low_High),dimnames=dimnames)
  RNA_seq_dataframe_Low_High=RNA_seq_dataframe_Low_High[rowMeans(RNA_seq_dataframe_Low_High)>1,]#过滤掉表达量低的基因

  conNum=ncol(RNA_seq_dataframe_Low)
  treatNum=ncol(RNA_seq_dataframe_High)
  Type_by_order=c(rep(1,conNum), rep(2,treatNum))

  #Differentially expressed genes analysis
  logFCfilter=0.585
  fdrFilter=0.05
  DEG_outTab_2vs3=data.frame()
  for(i in row.names(RNA_seq_dataframe_Low_High)){
    unigene_exp=data.frame(expression=RNA_seq_dataframe_Low_High[i,], Type=Type_by_order)
    wilcoxTest=wilcox.test(expression ~ Type, data=unigene_exp)
    pvalue=wilcoxTest$p.value
    conGeneMeans=mean(RNA_seq_dataframe_Low_High[i,1:conNum])
    treatGeneMeans=mean(RNA_seq_dataframe_Low_High[i,(conNum+1):ncol(RNA_seq_dataframe_Low_High)])
    logFC=treatGeneMeans-conGeneMeans
    conMed=median(RNA_seq_dataframe_Low_High[i,1:conNum])
    treatMed=median(RNA_seq_dataframe_Low_High[i,(conNum+1):ncol(RNA_seq_dataframe_Low_High)])
    diffMed=treatMed-conMed
    if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
      DEG_outTab_2vs3=rbind(DEG_outTab_2vs3,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
    }
  }
  pValue=DEG_outTab_2vs3[,"pValue"]
  fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
  DEG_outTab_2vs3=cbind(DEG_outTab_2vs3, fdr=fdr)
  Sig_DEG_outTab_2vs3=DEG_outTab_2vs3[( abs(as.numeric(as.vector(DEG_outTab_2vs3$logFC)))>logFCfilter & as.numeric(as.vector(DEG_outTab_2vs3$fdr))<fdrFilter),]

  #Getting intersection of DEGs and output the combined data set
  Sig_DEG_outTab_intersection <- intersect(Sig_DEG_outTab_1vs3$gene,Sig_DEG_outTab_2vs3$gene)
  inter_DEG_RNA_seq_dataframe=RNA_seq_dataframe[Sig_DEG_outTab_intersection,]
  inter_DEG_RNA_seq_dataframe=t(inter_DEG_RNA_seq_dataframe)
  sameSample=intersect(row.names(inter_DEG_RNA_seq_dataframe),row.names(MCKV_group))
  inter_DEG_RNA_seq_dataframe=inter_DEG_RNA_seq_dataframe[sameSample,]
  inter_DEG_clin_info=clin_info[sameSample,]
  inter_DEG_clin_info$futime=inter_DEG_clin_info$futime/365
  inter_DEG_clin_info_RNA_seq_dataframe=cbind(inter_DEG_clin_info,inter_DEG_RNA_seq_dataframe)

#Identifying prognostic related DEGs
coxPfilter=0.05
sel_DEG_outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(inter_DEG_clin_info_RNA_seq_dataframe[,3:ncol(inter_DEG_clin_info_RNA_seq_dataframe)])){
  #cox analysis
  cox <- coxph(Surv(futime, fustat) ~ inter_DEG_clin_info_RNA_seq_dataframe[,i], data = inter_DEG_clin_info_RNA_seq_dataframe)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  sigGenes=c(sigGenes,i)
  sel_DEG_outTab=rbind(sel_DEG_outTab,
                       cbind(id=i,
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR.95L=coxSummary$conf.int[,"lower .95"],
                             HR.95H=coxSummary$conf.int[,"upper .95"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}

sel_DEG_outTab$FDR <- p.adjust(sel_DEG_outTab$pvalue, method = "fdr")
sel_DEG_outTab_Filter <- sel_DEG_outTab[sel_DEG_outTab$FDR<0.25,]
}
