# MC-KV
We propose a Monte-Carlo K-Nearest Neighbor Voting (MC-KV) approach, a novel prognosis-oriented classifier based on semi-supervised learning for molecular subtyping of colorectal cancer. Based on the subtyping results of MC-KV method, we then establish an analysis platform for prognostic-related DEGs selection and prediction model construction.
## 1. Installation
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
 ## 2. Data preparation and input parameters
 * RNA_seq_dataframe: Data frame of the RNA sequencing profiles from public database<br>
 * expr_dataframe: Data frame of the biological input data (e.g. results from ssGSEA)<br>
 * clin_info: Data frame of the clinical information containing full survival status and survival time<br>
 * MCKV_group: the manually grouping results based on MCKNNV_tSNE<br>
 * candidate_gene: Data frame of the clinical information and RNA sequencing profiles of the candidate genes selected by WGCNA (top 5 prognostic DEGs in each module)<br>
 * clin_patho_parameter: Data frame of the clinicopathological parameters which are included in nomogram<br>
 * parameters:<br>
 p: the proportion of training set to total samples<br>
 n: the replicated times of kNN algorithm in the MC simulation<br>
 y: the cutoff year of manual labeling(e.g. 1, 3, 5)<br>
     y1: the first cutoff year of manual labeling<br>
     y2: the second cutoff year of manual labeling<br>
     y3: the third cutoff year of manual labeling<br>
 
 
 ## 3. Main procedure of MC-KV project
 **Step 1: Creating MC-KV algorithm based on survival status**<br> 
 MCKNNV(expr_dataframe, clin_info, p, n, y)<br> 
 **Step 2: Creating tSNE distribution based on MC-KV results**<br>
 MCKNNV_tSNE(expr_dataframe, clin_info, p, n, y1, y2, y3)<br> 
 **Step 3: Comparing survival difference between the groups identified by MC-KV**<br>
 Cluster_Surv_Diff(MCKV_group)<br>
 **Step 4: DEG analysis and genes intersection for subsequent process**<br>
 DEG_intersection(RNA_seq_dataframe, clin_info, MCKV_group)<br>
 **Step 5: Prognostic DEGs selection based on survival data**<br>
 prog_DEG_selection(RNA_seq_dataframe, clin_info, MCKV_group)<br>
 **Step 6: Run WGCNA with all intersected DEGs, and select top 5 prognostic DEGs in each module as the candidate genes**<br>
 **Step 7: Modeling gene selection and risk model construction based on DEGs from different groups**<br>
 Mod_Diff(candidate_gene)<br>
 **Step 8: Establishing nomogram based on the risk model and clinicopathological parameters**<br>
 Nomogram(candidate_gene,clin_patho_parameter)<br>
 
