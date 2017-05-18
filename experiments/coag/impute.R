library(dplyr)
library(mice)

#check missingness
coag %>% is.na %>% colSums

#create complete case analysis. almost halves data
#Coag <- AcitFactors %>% select(-hr0_apc) %>% na.omit

#use mice to impute
ini <- mice(coag, maxit = 0)
PredMatrix <- ini$predictorMatrix
PredMatrix[, c("acitnum")] <- 0
method <- ini$method
method[c("age", "bmi", "mechtype", "blunt", "tbi", "edarrivalgcs",
         "numribfxs", "admitday_intubated", "iss", "aishead1",
         "aisface2", "aischest3", "aisabdomen4", "aisextremity5",
         "prehosp_crystalloids", "icu_0to6h_blood_units",
         "minutesToEdArrival", "minutesToBloodDraw", "daysToDeath", "died")] <- ""
method[method != ""] <- "cart"

CoagImputed <- mice(coag, method = method, pred = PredMatrix,
                    minbucket = 5, m = 5, maxit = 50,  seed = 1991)

CoagImputed1 <- mice::complete(CoagImputed) %>% tbl_df
CoagImputed2 <- mice::complete(CoagImputed, action = 2) %>% tbl_df
CoagImputed3 <- mice::complete(CoagImputed, action = 3) %>% tbl_df
CoagImputed4 <- mice::complete(CoagImputed, action = 4) %>% tbl_df
CoagImputed5 <- mice::complete(CoagImputed, action = 5) %>% tbl_df

tegProteinImputed1 <- CoagImputed1 %>%
  select(hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab, hr0_plts,
         
         sample0h_crt_r, sample0h_crt_k, sample0h_crt_ma, sample0h_crt_ly30)

tegProteinImputed2 <- CoagImputed2 %>%
  select(hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab, hr0_plts,
         
         sample0h_crt_r, sample0h_crt_k, sample0h_crt_ma, sample0h_crt_ly30)

tegProteinImputed3 <- CoagImputed3 %>%
  select(hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab, hr0_plts,
         
         sample0h_crt_r, sample0h_crt_k, sample0h_crt_ma, sample0h_crt_ly30)

tegProteinImputed4 <- CoagImputed4 %>%
  select(hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab, hr0_plts,
         
         sample0h_crt_r, sample0h_crt_k, sample0h_crt_ma, sample0h_crt_ly30)

tegProteinImputed5 <- CoagImputed5 %>%
  select(hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab, hr0_plts,
         
         sample0h_crt_r, sample0h_crt_k, sample0h_crt_ma, sample0h_crt_ly30)

sparseCoagImputed1 <- CoagImputed1 %>%
  select(hr0_temp, hr0_hr, hr0_resprate, hr0_sbp,
         
         hr0_paco2, hr0_pao2, hr0_hco3, hr0_serumco2, hr0_bun, hr0_basedefexc, 
         
         hr0_wbc, hr0_hct, hr0_hgb, hr0_plts,
         
         hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab)
