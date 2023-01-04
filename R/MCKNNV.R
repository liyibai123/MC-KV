

#' @title MCKNNV
#' @description Creating MC-KNNV algorithm based on survival status
#' @details Input takes two data.frame expr_dataframe and clin_info with three variables p, n, and y
#' @param expr_dataframe Data frame of the biological input data (e.g. results from ssGSEA)
#' @param clin_info Data frame of the clinical information containing full survival status and survival time.
#' @param p the proportion of training set to total samples
#' @param n the replicated times of kNN algorithm in the MC simulation
#' @param y the cutoff year of manual labeling
#' @export


MCKNNV <- function(expr_dataframe, clin_info, p, n, y){
  # Manually labeling based on 3-year survival status
  clin_info$futime=clin_info$futime/365
  clin_info_yy <- clin_info
  clin_info_yy$LABEL[clin_info$fustat==1&clin_info$futime<=y] <- 1
  clin_info_yy$LABEL[clin_info$futime>y] <- 2
  clin_info_yy$LABEL[clin_info$fustat==0&clin_info$futime<=y] <- NA

  # Combining the labeled clinical information and expression data
  clin_info_yy_labeled <- clin_info_yy[!is.na(clin_info_yy$LABEL),]
  labeled_sample_yy=intersect(row.names(expr_dataframe),row.names(clin_info_yy_labeled))
  expr_dataframe_yy_labeled=expr_dataframe[labeled_sample_yy,]
  clin_info_yy_labeled <- clin_info_yy_labeled[labeled_sample_yy,]
  expr_dataframe_yy_labeled= cbind(clin_info_yy_labeled$LABEL,expr_dataframe_yy_labeled)
  colnames(expr_dataframe_yy_labeled)[1] <- "LABEL"

  # Combining the unlabeled clinical information and expression data
  clin_info_yy_unlabeled <- clin_info_yy[is.na(clin_info_yy$LABEL),]
  unlabeled_sample_yy=intersect(row.names(expr_dataframe),row.names(clin_info_yy_unlabeled))
  expr_dataframe_yy_unlabeled=expr_dataframe[unlabeled_sample_yy,]
  clin_info_yy_unlabeled <- clin_info_yy_unlabeled[unlabeled_sample_yy,]
  expr_dataframe_yy_unlabeled= cbind(clin_info_yy_unlabeled$LABEL,expr_dataframe_yy_unlabeled)
  colnames(expr_dataframe_yy_unlabeled)[1] <- "LABEL"

  # Running MC_KNNV process
  library(kknn)
  expr_dataframe_yy_labeled$LABEL <- factor(expr_dataframe_yy_labeled$LABEL,levels = c("1","2"),labels = c("1","2"))
  MC_KNNV_yy_result=data.frame(id=row.names(expr_dataframe))
  for(i in 1:n){
    print(i)
    index<-sample(1:nrow(expr_dataframe_yy_labeled),round(p*nrow(expr_dataframe_yy_labeled)))
    train<-expr_dataframe_yy_labeled[index,]
    test<-expr_dataframe_yy_labeled[-index,]
    test <- rbind(test,expr_dataframe_yy_unlabeled)
    pre<-kknn(LABEL~.,train,test, k=7, distance = 2, kernel= c("triangular"))
    pred<-data.frame(id=row.names(test), pred=fitted(pre))
    colnames(pred)[2]=paste0("pred",i)
    MC_KNNV_yy_result=merge(MC_KNNV_yy_result,pred,by="id",all = T)
  }
    rownames(MC_KNNV_yy_result)=MC_KNNV_yy_result$id
    MC_KNNV_yy_result=t(MC_KNNV_yy_result)[-1,]

    vot_yy <- data.frame(sampID=colnames(MC_KNNV_yy_result), Votdata_yy=NA)

    for(i in 1:ncol(MC_KNNV_yy_result)){
      voting_yy <- as.data.frame(table(MC_KNNV_yy_result[,i]))
      voting_yy$Var1 <- as.numeric(as.character(voting_yy$Var1))
      if(nrow(voting_yy)<2){
        added = data.frame(Var1=3-voting_yy$Var1,Freq=0)
        voting_yy <- rbind(voting_yy, added)
      }
      G1 <- voting_yy$Freq[voting_yy$Var1==1]
      G2 <- voting_yy$Freq[voting_yy$Var1==2]
      X <- (G1-G2)/(G1+G2)
      Y <- G1/(G1+G2)
      vot_yy$Votdata_yy[i] <-X
      vot_yy$Votrate_yy[i] <-Y
    }

    rownames(vot_yy) <- vot_yy$sampID
    vot_yy$cat <- cut(vot_yy$Votdata_yy, c(-10,0,10))

    comrowname <- intersect(rownames(vot_yy), rownames(clin_info))
    tab1<- vot_yy[comrowname,]
    tab2<- clin_info[comrowname,]

    survtab_yy <- merge(tab1, tab2, by= "row.names")

  }


