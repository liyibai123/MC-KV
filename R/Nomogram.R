#' @title Nomogram
#' @description Establishing nomogram based on the gene model and clinicopathological parameters
#' @details Input takes data frame of candidate genes created by the exports of DEG_selection and WGCNA, with the clinicopathological parameters downloaded from public database
#' @param candidate_gene Data frame of the clinical information and RNA sequencing profiles of the candidate genes selected by WGCNA (top 5 prognostic DEGs in each module)
#' @param clin_patho_parameter Data frame of the clinicopathological parameters which are included in nomogram
#' @export


Nomogram <- function(candidate_gene,clin_patho_parameter){
  library("VRPM")
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

trainRisk=trainRisk[,c("futime", "fustat", "risk")]
clin_patho_parameter=clin_patho_parameter[apply(clin_patho_parameter,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
clin_patho_parameter$Age=as.numeric(clin_patho_parameter$Age)

#Combining gene model and clinicopathological parameters
samSample=intersect(row.names(trainRisk), row.names(clin_patho_parameter))
trainRisk=trainRisk[samSample,,drop=F]
clin_patho_parameter=clin_patho_parameter[samSample,,drop=F]
nomo_data=cbind(trainRisk, clin_patho_parameter)

#Establishing nomogram
library(regplot)
res.cox=coxph(Surv(futime, fustat) ~ risk + Age + Stage, data = nomo_data)
output = regplot(res.cox,
        plots = c("density", "boxes"),
        clickable=F,
        title="",
        points=TRUE,
        droplines=TRUE,
        observation=nomo_data[1,],
        rank="sd",
        failtime = c(1,3,5),
        prfail = T)
}
