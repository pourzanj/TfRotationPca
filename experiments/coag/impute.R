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
                    minbucket = 5, m = 1, maxit = 50,  seed = 1991)

CoagImputed1 <- mice::complete(CoagImputed) %>% tbl_df
