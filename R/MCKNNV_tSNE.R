
#' @title MCKNNV_tSNE
#' @description Creating tSNE distribution based on the MC-KNNV results
#' @details Input takes several data frames which are the exports of MCKNNV with different y
#' @param expr_dataframe Data frame of the biological input data (e.g. results from ssGSEA)
#' @param clin_info Data frame of the clinical information containing full survival status and survival time.
#' @param p the proportion of training set to total samples
#' @param n the replicated times of kNN algorithm in the MC simulation
#' @param y1 the first cutoff year of manual labeling
#' @param y2 the second cutoff year of manual labeling
#' @param y3 the third cutoff year of manual labeling
#' @export


MCKNNV_tSNE <- function(expr_dataframe, clin_info, p, n, y1, y2, y3){
#output for y1
  # Manually labeling based on survival status
  clin_info$futime=clin_info$futime/365
  clin_info_y1y <- clin_info
  clin_info_y1y$LABEL[clin_info$fustat==1&clin_info$futime<=y1] <- 1
  clin_info_y1y$LABEL[clin_info$futime>y1] <- 2
  clin_info_y1y$LABEL[clin_info$fustat==0&clin_info$futime<=y1] <- NA

  # Combining the labeled clinical information and expression data
  clin_info_y1y_labeled <- clin_info_y1y[!is.na(clin_info_y1y$LABEL),]
  labeled_sample_y1y=intersect(row.names(expr_dataframe),row.names(clin_info_y1y_labeled))
  expr_dataframe_y1y_labeled=expr_dataframe[labeled_sample_y1y,]
  clin_info_y1y_labeled <- clin_info_y1y_labeled[labeled_sample_y1y,]
  expr_dataframe_y1y_labeled= cbind(clin_info_y1y_labeled$LABEL,expr_dataframe_y1y_labeled)
  colnames(expr_dataframe_y1y_labeled)[1] <- "LABEL"

  # Combining the unlabeled clinical information and expression data
  clin_info_y1y_unlabeled <- clin_info_y1y[is.na(clin_info_y1y$LABEL),]
  unlabeled_sample_y1y=intersect(row.names(expr_dataframe),row.names(clin_info_y1y_unlabeled))
  expr_dataframe_y1y_unlabeled=expr_dataframe[unlabeled_sample_y1y,]
  clin_info_y1y_unlabeled <- clin_info_y1y_unlabeled[unlabeled_sample_y1y,]
  expr_dataframe_y1y_unlabeled= cbind(clin_info_y1y_unlabeled$LABEL,expr_dataframe_y1y_unlabeled)
  colnames(expr_dataframe_y1y_unlabeled)[1] <- "LABEL"

  # Running MC_KNNV process
  library(kknn)
  expr_dataframe_y1y_labeled$LABEL <- factor(expr_dataframe_y1y_labeled$LABEL,levels = c("1","2"),labels = c("1","2"))
  MC_KNNV_y1y_result=data.frame(id=row.names(expr_dataframe))
  for(i in 1:n){
    print(i)
    index<-sample(1:nrow(expr_dataframe_y1y_labeled),round(p*nrow(expr_dataframe_y1y_labeled)))
    train<-expr_dataframe_y1y_labeled[index,]
    test<-expr_dataframe_y1y_labeled[-index,]
    test <- rbind(test,expr_dataframe_y1y_unlabeled)
    pre<-kknn(LABEL~.,train,test, k=7, distance = 2, kernel= c("triangular"))
    pred<-data.frame(id=row.names(test), pred=fitted(pre))
    colnames(pred)[2]=paste0("pred",i)
    MC_KNNV_y1y_result=merge(MC_KNNV_y1y_result,pred,by="id",all = T)
  }
  rownames(MC_KNNV_y1y_result)=MC_KNNV_y1y_result$id
  MC_KNNV_y1y_result=t(MC_KNNV_y1y_result)[-1,]

  vot_y1y <- data.frame(sampID=colnames(MC_KNNV_y1y_result), Votdata_y1y=NA)

  for(i in 1:ncol(MC_KNNV_y1y_result)){
    voting_y1y <- as.data.frame(table(MC_KNNV_y1y_result[,i]))
    voting_y1y$Var1 <- as.numeric(as.character(voting_y1y$Var1))
    if(nrow(voting_y1y)<2){
      added = data.frame(Var1=3-voting_y1y$Var1,Freq=0)
      voting_y1y <- rbind(voting_y1y, added)
    }
    G1 <- voting_y1y$Freq[voting_y1y$Var1==1]
    G2 <- voting_y1y$Freq[voting_y1y$Var1==2]
    X <- (G1-G2)/(G1+G2)
    Y <- G1/(G1+G2)
    vot_y1y$Votdata_y1y[i] <-X
    vot_y1y$Votrate_y1y[i] <-Y
  }

  rownames(vot_y1y) <- vot_y1y$sampID
  vot_y1y$cat <- cut(vot_y1y$Votdata_y1y, c(-10,0,10))

  comrowname <- intersect(rownames(vot_y1y), rownames(clin_info))
  tab1<- vot_y1y[comrowname,]
  tab2<- clin_info[comrowname,]

  survtab_y1y <- merge(tab1, tab2, by= "row.names")

#output for y2
  # Manually labeling based on survival status
  clin_info_y2y <- clin_info
  clin_info_y2y$LABEL[clin_info$fustat==1&clin_info$futime<=y2] <- 1
  clin_info_y2y$LABEL[clin_info$futime>y2] <- 2
  clin_info_y2y$LABEL[clin_info$fustat==0&clin_info$futime<=y2] <- NA

  # Combining the labeled clinical information and expression data
  clin_info_y2y_labeled <- clin_info_y2y[!is.na(clin_info_y2y$LABEL),]
  labeled_sample_y2y=intersect(row.names(expr_dataframe),row.names(clin_info_y2y_labeled))
  expr_dataframe_y2y_labeled=expr_dataframe[labeled_sample_y2y,]
  clin_info_y2y_labeled <- clin_info_y2y_labeled[labeled_sample_y2y,]
  expr_dataframe_y2y_labeled= cbind(clin_info_y2y_labeled$LABEL,expr_dataframe_y2y_labeled)
  colnames(expr_dataframe_y2y_labeled)[1] <- "LABEL"

  # Combining the unlabeled clinical information and expression data
  clin_info_y2y_unlabeled <- clin_info_y2y[is.na(clin_info_y2y$LABEL),]
  unlabeled_sample_y2y=intersect(row.names(expr_dataframe),row.names(clin_info_y2y_unlabeled))
  expr_dataframe_y2y_unlabeled=expr_dataframe[unlabeled_sample_y2y,]
  clin_info_y2y_unlabeled <- clin_info_y2y_unlabeled[unlabeled_sample_y2y,]
  expr_dataframe_y2y_unlabeled= cbind(clin_info_y2y_unlabeled$LABEL,expr_dataframe_y2y_unlabeled)
  colnames(expr_dataframe_y2y_unlabeled)[1] <- "LABEL"

  # Running MC_KNNV process
  library(kknn)
  expr_dataframe_y2y_labeled$LABEL <- factor(expr_dataframe_y2y_labeled$LABEL,levels = c("1","2"),labels = c("1","2"))
  MC_KNNV_y2y_result=data.frame(id=row.names(expr_dataframe))
  for(i in 1:n){
    print(i)
    index<-sample(1:nrow(expr_dataframe_y2y_labeled),round(p*nrow(expr_dataframe_y2y_labeled)))
    train<-expr_dataframe_y2y_labeled[index,]
    test<-expr_dataframe_y2y_labeled[-index,]
    test <- rbind(test,expr_dataframe_y2y_unlabeled)
    pre<-kknn(LABEL~.,train,test, k=7, distance = 2, kernel= c("triangular"))
    pred<-data.frame(id=row.names(test), pred=fitted(pre))
    colnames(pred)[2]=paste0("pred",i)
    MC_KNNV_y2y_result=merge(MC_KNNV_y2y_result,pred,by="id",all = T)
  }
  rownames(MC_KNNV_y2y_result)=MC_KNNV_y2y_result$id
  MC_KNNV_y2y_result=t(MC_KNNV_y2y_result)[-1,]

  vot_y2y <- data.frame(sampID=colnames(MC_KNNV_y2y_result), Votdata_y2y=NA)

  for(i in 1:ncol(MC_KNNV_y2y_result)){
    voting_y2y <- as.data.frame(table(MC_KNNV_y2y_result[,i]))
    voting_y2y$Var1 <- as.numeric(as.character(voting_y2y$Var1))
    if(nrow(voting_y2y)<2){
      added = data.frame(Var1=3-voting_y2y$Var1,Freq=0)
      voting_y2y <- rbind(voting_y2y, added)
    }
    G1 <- voting_y2y$Freq[voting_y2y$Var1==1]
    G2 <- voting_y2y$Freq[voting_y2y$Var1==2]
    X <- (G1-G2)/(G1+G2)
    Y <- G1/(G1+G2)
    vot_y2y$Votdata_y2y[i] <-X
    vot_y2y$Votrate_y2y[i] <-Y
  }

  rownames(vot_y2y) <- vot_y2y$sampID
  vot_y2y$cat <- cut(vot_y2y$Votdata_y2y, c(-10,0,10))

  comrowname <- intersect(rownames(vot_y2y), rownames(clin_info))
  tab1<- vot_y2y[comrowname,]
  tab2<- clin_info[comrowname,]

  survtab_y2y <- merge(tab1, tab2, by= "row.names")

#output for y3
  # Manually labeling based on survival status
  clin_info_y3y <- clin_info
  clin_info_y3y$LABEL[clin_info$fustat==1&clin_info$futime<=y3] <- 1
  clin_info_y3y$LABEL[clin_info$futime>y3] <- 2
  clin_info_y3y$LABEL[clin_info$fustat==0&clin_info$futime<=y3] <- NA

  # Combining the labeled clinical information and expression data
  clin_info_y3y_labeled <- clin_info_y3y[!is.na(clin_info_y3y$LABEL),]
  labeled_sample_y3y=intersect(row.names(expr_dataframe),row.names(clin_info_y3y_labeled))
  expr_dataframe_y3y_labeled=expr_dataframe[labeled_sample_y3y,]
  clin_info_y3y_labeled <- clin_info_y3y_labeled[labeled_sample_y3y,]
  expr_dataframe_y3y_labeled= cbind(clin_info_y3y_labeled$LABEL,expr_dataframe_y3y_labeled)
  colnames(expr_dataframe_y3y_labeled)[1] <- "LABEL"

  # Combining the unlabeled clinical information and expression data
  clin_info_y3y_unlabeled <- clin_info_y3y[is.na(clin_info_y3y$LABEL),]
  unlabeled_sample_y3y=intersect(row.names(expr_dataframe),row.names(clin_info_y3y_unlabeled))
  expr_dataframe_y3y_unlabeled=expr_dataframe[unlabeled_sample_y3y,]
  clin_info_y3y_unlabeled <- clin_info_y3y_unlabeled[unlabeled_sample_y3y,]
  expr_dataframe_y3y_unlabeled= cbind(clin_info_y3y_unlabeled$LABEL,expr_dataframe_y3y_unlabeled)
  colnames(expr_dataframe_y3y_unlabeled)[1] <- "LABEL"

  # Running MC_KNNV process
  library(kknn)
  expr_dataframe_y3y_labeled$LABEL <- factor(expr_dataframe_y3y_labeled$LABEL,levels = c("1","2"),labels = c("1","2"))
  MC_KNNV_y3y_result=data.frame(id=row.names(expr_dataframe))
  for(i in 1:n){
    print(i)
    index<-sample(1:nrow(expr_dataframe_y3y_labeled),round(p*nrow(expr_dataframe_y3y_labeled)))
    train<-expr_dataframe_y3y_labeled[index,]
    test<-expr_dataframe_y3y_labeled[-index,]
    test <- rbind(test,expr_dataframe_y3y_unlabeled)
    pre<-kknn(LABEL~.,train,test, k=7, distance = 2, kernel= c("triangular"))
    pred<-data.frame(id=row.names(test), pred=fitted(pre))
    colnames(pred)[2]=paste0("pred",i)
    MC_KNNV_y3y_result=merge(MC_KNNV_y3y_result,pred,by="id",all = T)
  }
  rownames(MC_KNNV_y3y_result)=MC_KNNV_y3y_result$id
  MC_KNNV_y3y_result=t(MC_KNNV_y3y_result)[-1,]

  vot_y3y <- data.frame(sampID=colnames(MC_KNNV_y3y_result), Votdata_y3y=NA)

  for(i in 1:ncol(MC_KNNV_y3y_result)){
    voting_y3y <- as.data.frame(table(MC_KNNV_y3y_result[,i]))
    voting_y3y$Var1 <- as.numeric(as.character(voting_y3y$Var1))
    if(nrow(voting_y3y)<2){
      added = data.frame(Var1=3-voting_y3y$Var1,Freq=0)
      voting_y3y <- rbind(voting_y3y, added)
    }
    G1 <- voting_y3y$Freq[voting_y3y$Var1==1]
    G2 <- voting_y3y$Freq[voting_y3y$Var1==2]
    X <- (G1-G2)/(G1+G2)
    Y <- G1/(G1+G2)
    vot_y3y$Votdata_y3y[i] <-X
    vot_y3y$Votrate_y3y[i] <-Y
  }

  rownames(vot_y3y) <- vot_y3y$sampID
  vot_y3y$cat <- cut(vot_y3y$Votdata_y3y, c(-10,0,10))

  comrowname <- intersect(rownames(vot_y3y), rownames(clin_info))
  tab1<- vot_y3y[comrowname,]
  tab2<- clin_info[comrowname,]

  survtab_y3y <- merge(tab1, tab2, by= "row.names")


#combining the results of y1, y2 and y3
  survtab_merge <- merge(survtab_y1y,survtab_y2y,by="sampID",suffixes = c(".y1y",".y2y"),sort=T)
  survtab_merge <- merge(survtab_merge,survtab_y3y,by="sampID",sort=T)
  library(tidyverse)
  Votdata_merge <- survtab_merge %>% select(1,3,9,15)

  library(Rtsne)
  rownames(Votdata_merge) <- Votdata_merge$sampID
  Votdata_merge1=Votdata_merge[,-1]

  data_tsne = Rtsne(Votdata_merge, dims=2, perplexity=10, max_iter=2000)
  plotdata <- as.data.frame(data_tsne$Y)
  plotdata_manu <- plotdata

  pic_tsne <- ggplot(plotdata_manu,
                     aes(x=V1, y=V2))+
    geom_point(size=2.5)  +
    theme_bw()
}
