# Empirical Likelihood Inference for Area under the ROC Curve using Ranked Set Samples

The area under a receiver operating characteristic curve (AUC) is a useful tool to assess the performance of continuous-scale diagnostic tests on binary classification. 
In [Moon et al.](https://arxiv.org/abs/2010.12185), we propose an empirical likelihood (EL) method to construct confidence intervals for the AUC from data 
collected by ranked set sampling (RSS). The proposed EL-based method enables inferences without assumptions required in existing nonparametric methods and 
takes advantage of the sampling efficiency of RSS. We show that for both balanced and unbalanced RSS, the EL-based point estimate is the Mann-Whitney statistic, 
and confidence intervals can be obtained from a scaled chi-square distribution. 

Simulation studies and two case studies on diabetes and chronic kidney disease data are carried and show that the proposed method performs better than the existing methods.

# Data Availability
Two case studies are based upon data collected from [National Health and Nutrition Examination Survey (NHANES)](https://www.cdc.gov/nchs/nhanes/). 
The NHANES data are publicly available and retrieved in R.

# Reproducing Simulation and Case Studies
1. Simulation studies \
Run `run_simulation_BRSS.R` and `run_simulation_URSS.R` for balanced RSS and unbalanced RSS studies. For plots, run `results_BRSS.R` and `results_URSS.R`. Run `run_simulation_BRSS_reference.R` for comparing reference distributions.
2. Diabetes case study \
Run `run_diabetes_application.R` and `results_diabetes.R` for applications and plots, respectively.
3. Chronic kidney disease case study \
Run `run_sim_ckd_application.R` and `results_ckd.R` for applications and plots, respectively. To retrieve the CKD data used in this study `data_ckd.Rdata`, run `data_ckd.R`.
4. Simulation time \
Run `simulation_time.R` for computing the simulation time.
6. Toy example \
Run `example.R` for a toy example that computes confidence intervals of BRSS-EL and URSS-EL.
