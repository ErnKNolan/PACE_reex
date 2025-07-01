#AUTHOR: Erin Nolan
#DATE CREATED: 18 FEB 2025
#PURPOSE: Import and data management

pacman::p_load(tidyverse,here)

pace_noninf <- read.csv(here("Data","Original","pace_wide_noninf.csv"))
pace_noninf <- pace_noninf %>%
  mutate(LD = factor(dose_fin),
         ses = ifelse(Decile < 6, "Most disadvantaged","Least disadvantaged"),
         RA_name = ifelse(RA_name %in% c("Major Cities of Australia"),"Major Cities","Regional/Remote"),
         stype = factor(stype))

#The second PACE trial
pace_sup <- read.csv(here("Data","Original","PACETotalPAT1T2_RAW.csv"))
pace_sup2 <- pace_sup %>%
  #name variables and groupings the same as the 3rd pace trial
  mutate(LD = factor(ifelse(intcont == 2,"Ctrl","H")),
         ses = factor(ifelse(SEIFACAT16 == "2. Least Disadvantaged [NSW 2016]","Least disadvantaged","Most disadvantaged")),
         RA_name = factor(ifelse(RemoteDICHOM16 == "2. Inner/Outer Regional/Remote Australia","Regional/Remote","Major Cities")),
         STYPE = factor(STYPE)) %>%
  rename(tname = teacher_IDFIN,
         schoolid = SCHOOLID,
         stype = STYPE) %>%
  group_by(tname) %>%
  arrange(tname,time) %>%
  #keep 1 row per teacher, to later merge in wide format
  filter(row_number()==1) %>%
  dplyr::select(tname,schoolid,stype,Allocation_Date,LD,ses,RA_name)
  
#make a wide format of the two times
pace_wide <- pace_sup %>%
  pivot_wider(id_cols = teacher_IDFIN, names_from=time, values_from=total_paCTA) %>%
  rename(base_total_pa_fin = `1`,
         fu_total_pa_fin = `2`)

#combine the wide format outcomes with the rest of the data
pace_sup3 <- left_join(pace_sup2,pace_wide,by=c("tname" = "teacher_IDFIN"))


#merge the two datasets together
#the 2nd and 3rd pace trials already have the school id and teacher id spaced out, no need to redo IDs
#schoolid is in order of allocation date
pace <- bind_rows(pace_sup3,pace_noninf) %>%
  dplyr::select(-dose_fin,-Decile,-Allocation_Date) %>%
  mutate(schoolid = factor(schoolid),
         tname = factor(tname),
         RA_name = factor(RA_name),
         ses = factor(ses),
         trial = ifelse(as.numeric(as.character(schoolid)) <= 62,1,0))

saveRDS(pace,here("Data","pace.RDS"))

