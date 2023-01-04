#' @title Cluster_Surv_Diff
#' @description Comparing the survival difference between the group identified by MC-KNNV
#' @details Input takes the output of MCKNNV_tSNE
#' @param MCKV_group the manually grouping results based on MCKNNV_tSNE
#' @export


Cluster_Surv_Diff <- function(MCKV_group){
  library(pROC)
  library(survival)
  library(survminer)
  diff=survdiff(Surv(futime, fustat) ~cluster, data = MCKV_group)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ cluster, data = MCKV_group)

  surPlot=ggsurvplot(fit,
                     data=MCKV_group,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="cluster",
                     legend.labs=c("1","2","3"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("#B91425","#FEE407","#244999"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
}
