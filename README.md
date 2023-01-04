# MC-KV
We propose a Monte-Carlo K-Nearest Neighbor Voting (MC-KV) approach, a novel prognosis-oriented classifier based on semi-supervised learning for molecular subtyping of colorectal cancer. Based on the subtyping results of MC-KV method, we then establish an analysis platform for prognostic-related DEGs selection and prediction model construction.
# 1. Installation
Before using the MC-KV classifier and subsequent analysis methods, several dependent packages need to be downloaded and installed in advance: 
 ```install.package(kknn)   
 install.package(tidyverse)  
 install.package(Rtsne)   
 install.package(pROC)    
 install.package(survival)    
 install.package(survminer)   
 install.package(glmnet)    
 install.package(randomForestSRC)   
 install.package(VRPM)    
 install.package(regplot)
 ```
 
 # 2. Main procedure of MC-KV project
 **Step 1**: MCKNNV(expr_dataframe, clin_info, p, n, y)<br> 
 **Step 2**: MCKNNV_tSNE(expr_dataframe, clin_info, p, n, y1, y2, y3)<br> 
 **Step 3**: Cluster_Surv_Diff(MCKV_group)<br>
 **Step 4**: DEG_intersection(RNA_seq_dataframe, clin_info, MCKV_group)<br>
 **Step 5**: prog_DEG_selection(RNA_seq_dataframe, clin_info, MCKV_group)<br>
 **Step 6**: run WGCNA with all intersected DEGs, and select the top 5 prognostic DEGs in each module as the candidate genes<br>
 **Step 7**: Mod_Diff(candidate_gene)<br>
 **Step 8**: Nomogram(candidate_gene,clin_patho_parameter)<br>
 
 
