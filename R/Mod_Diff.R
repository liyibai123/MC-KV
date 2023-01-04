#' @title Mod_Diff
#' @description Modeling gene selection and model construction based on the DEGs from different groups
#' @details Input takes data frame of candidate genes created by the exports of DEG_selection and WGCNA
#' @param candidate_gene Data frame of the clinical information and RNA sequencing profiles of the candidate genes selected by WGCNA (top 5 prognostic DEGs in each module)
#' @export


Mod_Diff <- function(candidate_gene){
  library("glmnet")
  library("survival")
  candidate_gene$futime[candidate_gene$futime<=0]=0.003

  #Screening model genes from candidate genes by LASSO
  x=as.matrix(candidate_gene[,c(3:ncol(candidate_gene))])
  y=data.matrix(Surv(candidate_gene$futime, candidate_gene$fustat))
  fit=glmnet(x, y, family="cox", maxit=1000)
  cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
  coef=coef(fit, s=0.05)
  index=which(coef != 0)
  actCoef=coef[index]
  lassoGene=row.names(coef)[index]
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)

  #Constructing Random Survival Forrest model
  library(randomForestSRC)
  model_gene=candidate_gene[,c("futime","fustat",lassoGene)]
  rfsrc_pbcmy <- rfsrc(Surv(futime, fustat) ~ .,data = model_gene, nsplit = 10,  na.action = "na.impute", tree.err = TRUE,importance = TRUE)
  probility_rfsrc_pbcmy <- predict(rfsrc_pbcmy,model_gene)
  trainScore=probility_rfsrc_pbcmy$predicted
  outCol=c("futime","fustat",lassoGene)
  risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
  trainRisk=cbind(model_gene[,outCol],riskScore=as.vector(trainScore),risk)
}



