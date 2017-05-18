library(dplyr)
library(readr)
library(lubridate)

coag <- read_csv("experiments/coag/ACITdataset_UCSB_23Feb15.csv") %>%
  tbl_df() %>%
  select(acitnum, male, age, bmi,
         mechtype, blunt, tbi, edarrivalgcs, numribfxs, admitday_intubated,
         iss, aishead1, aisface2, aischest3, aisabdomen4, aisextremity5,
         
         prehosp_crystalloids, icu_0to6h_blood_units,
         
         injurydatetime, edarrivaldatetime, hr0_datetime, datetimeofdeath,
         
         hr0_temp, hr0_hr, hr0_resprate, hr0_sbp,
         
        hr0_paco2, hr0_pao2, hr0_hco3, hr0_serumco2, hr0_bun, hr0_basedefexc, 
        
        hr0_wbc, hr0_hct, hr0_hgb, hr0_plts,
         
         hr0_factorii, hr0_factorv, hr0_factorvii, hr0_factorviii,
         hr0_factorix, hr0_factorx, hr0_ddimer, hr0_pc,
         hr0_atiii, hr0_fibrinogen_rlab,
         
        sample0h_crt_r, sample0h_crt_k, sample0h_crt_ma, sample0h_crt_ly30) %>%
        
         #hr0_pt, hr0_ptt) %>%
  
  #clean some of the columns
  mutate(mechtype = as.character(mechtype)) %>%
  mutate(mechtype = ifelse(mechtype == "", as.character(NA), mechtype)) %>%
  mutate(mechtype = as.factor(mechtype)) %>%
  
  mutate(blunt = as.character(blunt)) %>%
  mutate(blunt = ifelse(blunt == "", as.character(NA), blunt)) %>%
  mutate(blunt = as.factor(blunt)) %>%
  
  #mutate(hr0_lactate = ifelse(hr0_lactate == "", as.character(NA), hr0_lactate)) %>%
  #mutate(hr0_lactate = ifelse(hr0_lactate == ">15.0", NA, hr0_lactate)) %>%
  #mutate(hr0_lactate = as.numeric(hr0_lactate)) %>%
  
  #wrangle times
  mutate_each(funs(parse_date_time(., "m/d/y H:M")),
              injurydatetime, edarrivaldatetime, hr0_datetime, datetimeofdeath) %>%
  
  mutate(minutesToEdArrival = (injurydatetime %--% edarrivaldatetime) / dminutes(1)) %>%
  
  #fix cases where minutes to arrival doesn't make sense because its negative
  mutate(minutesToEdArrival = ifelse(minutesToEdArrival <= 0, NA, minutesToEdArrival)) %>%
  mutate(minutesToBloodDraw = (edarrivaldatetime %--% hr0_datetime) / dminutes(1)) %>%
  mutate(minutesToBloodDraw = ifelse(!between(minutesToBloodDraw, -50, 90), NA, minutesToBloodDraw)) %>%
  mutate(daysToDeath = (edarrivaldatetime %--% datetimeofdeath) / ddays(1)) %>%
  mutate(died = !is.na(datetimeofdeath)) %>%
  select(-injurydatetime, -hr0_datetime, -edarrivaldatetime, -datetimeofdeath) %>%
  
  filter(iss >= 25) %>%
  filter(age <= 45) %>%
  filter(mechtype != "Other") %>%
  filter(mechtype != "Found down")
  
  #filter(mechtype == "MVC")