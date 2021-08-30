library(tidyverse)

# functions for NHANES data retrieval
source("./functions_ckd_data.R")

nhanesdat = function(year,alphabet) {
	demo = nhanes_load_data(paste0("DEMO_",alphabet), year) # race and gender
	urine = nhanes_load_data(paste0("ALB_CR_",alphabet), year) # ratio
	biopro = nhanes_load_data(paste0("BIOPRO_",alphabet), year) # creatine serum
	diabetes = nhanes_load_data(paste0("DIQ_",alphabet), year) # diabetes
	
	dat = biopro %>% 
		left_join(demo,by="SEQN") %>%
		left_join(urine,by="SEQN") %>%
		left_join(diabetes,by="SEQN") %>%
		na_if(".") %>%
		select(SEQN,RIAGENDR,RIDAGEYR,RIDRETH1,LBXSCR,URDACT,DIQ010) %>%
		drop_na() %>%
		mutate(FEMALE=case_when(RIAGENDR==2 ~ 1.018,
														RIAGENDR==1 ~ 1),
					 BLACK=case_when(RIDRETH1==4 ~ 1.159,
					 								TRUE ~ 1),
					 FEMALEA=case_when(RIAGENDR==2 ~ -0.329,
					 									RIAGENDR==1 ~ -0.411),
					 FEMALEK=case_when(RIAGENDR==2 ~ 0.7,
					 									RIAGENDR==1 ~ 0.9),
					 DIABETES=case_when(DIQ010==1 ~ 1,
					 									 TRUE ~ 0)) %>%
		filter(DIABETES==1,RIDAGEYR>20,RIDAGEYR<80) %>%
		mutate(MINSCR = pmin(LBXSCR/FEMALEK,1), MAXSCR = pmax(LBXSCR/FEMALEK,1)) %>%
		mutate(GFR=141*MINSCR^(FEMALEA)*MAXSCR^(-1.209)*0.993^(RIDAGEYR)*FEMALE*BLACK) %>% # GFR
		mutate(CKD=case_when( (GFR<60 | URDACT>30)  ~ 1,TRUE ~ 0)) # CKD
	return(dat)
}

# retrieve NHANES data
dat2017=nhanesdat("2017-2018","J")
dat2015=nhanesdat("2015-2016","I")
dat2013=nhanesdat("2013-2014","H")
dat2011=nhanesdat("2011-2012","G")
dat2009=nhanesdat("2009-2010","F")

# aggregated data
data.ckd = rbind(dat2009,dat2011,dat2013,dat2015,dat2017)

save(data.ckd, file="./data_ckd.Rdata")
