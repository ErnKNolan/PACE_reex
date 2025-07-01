pacman::p_load(mice,brms,here,tidyverse,marginaleffects,blockrand)

'%!in%' <- function(x,y)!('%in%'(x,y))
set.seed(866380771)
pace_preds <- readRDS(here("Data","pace_preds.RDS"))
imp <- readRDS(here("Data","imp.RDS"))

#Interim 1
#The first third of schools
in2_time1_pre <- imp %>% filter(as.numeric(as.character(schoolid))<= 37) %>%
  mutate(schoolid2 = as.numeric(as.character(schoolid)))

#Randomisation scheduling
set.seed(866380771)
MajorCities_s1_rand <- blockrand(n=109,
                              num.levels=3,
                              levels = c("Ctrl","L","H"),
                              stratum = c("MajorCities.s1"),
                              block.sizes = c(1))%>%
  mutate(rownumber=row_number())
set.seed(866380771)
MajorCities_s2_rand <- blockrand(n=109,
                                 num.levels=3,
                                 levels = c("Ctrl","L","H"),
                                 stratum = c("MajorCities.s2"),
                                 block.sizes = c(1))%>%
  mutate(rownumber=row_number())
set.seed(866380771)
MajorCities_s3_rand <- blockrand(n=109,
                                 num.levels=3,
                                 levels = c("Ctrl","L","H"),
                                 stratum = c("MajorCities.s3"),
                                 block.sizes = c(1))%>%
  mutate(rownumber=row_number())
set.seed(866380771)
RegionalRemote_s1_rand <- blockrand(n=109,
                                 num.levels=3,
                                 levels = c("Ctrl","L","H"),
                                 stratum = c("Regional/Remote.s1"),
                                 block.sizes = c(1)) %>%
  mutate(rownumber=row_number())
set.seed(866380771)
RegionalRemote_s2_rand <- blockrand(n=109,
                                 num.levels=3,
                                 levels = c("Ctrl","L","H"),
                                 stratum = c("Regional/Remote.s2"),
                                 block.sizes = c(1)) %>%
  mutate(rownumber=row_number())
set.seed(866380771)
RegionalRemote_s3_rand <- blockrand(n=109,
                                 num.levels=3,
                                 levels = c("Ctrl","L","H"),
                                 stratum = c("Regional/Remote.s3"),
                                 block.sizes = c(1)) %>%
  mutate(rownumber=row_number())

#get the school ids and order them by randomisation order
schoolid_rand <- in2_time1_pre %>% 
  group_by(schoolid) %>% 
  filter(row_number()==1) %>% 
  ungroup()%>%
  dplyr::select(schoolid,RA_name,stype) %>%
  mutate(schoolid = as.numeric(as.character(schoolid))) %>%
  arrange(schoolid)

#assign the groups by stratification
#major cities
#s1
school_majorcities_s1_rand <- schoolid_rand %>% 
  filter(RA_name == "Major Cities",stype=="1") %>%
  mutate(rownumber = row_number()) %>%
  left_join(MajorCities_s1_rand,by="rownumber")
#s2
school_majorcities_s2_rand <- schoolid_rand %>% 
  filter(RA_name == "Major Cities",stype=="2") %>%
  mutate(rownumber = row_number()) %>%
  left_join(MajorCities_s2_rand,by="rownumber")
#s3
school_majorcities_s3_rand <- schoolid_rand %>% 
  filter(RA_name == "Major Cities",stype=="3") %>%
  mutate(rownumber = row_number()) %>%
  left_join(MajorCities_s3_rand,by="rownumber")

#regional remote
#s1
school_regremote_s1_rand <- schoolid_rand %>% 
  filter(RA_name == "Regional/Remote",stype=="1") %>%
  mutate(rownumber = row_number()) %>%
  left_join(RegionalRemote_s1_rand,by="rownumber")
#s2
school_regremote_s2_rand <- schoolid_rand %>% 
  filter(RA_name == "Regional/Remote",stype=="2") %>%
  mutate(rownumber = row_number()) %>%
  left_join(RegionalRemote_s2_rand,by="rownumber")
#s3
school_regremote_s3_rand <- schoolid_rand %>% 
  filter(RA_name == "Regional/Remote",stype=="3") %>%
  mutate(rownumber = row_number()) %>%
  left_join(RegionalRemote_s3_rand,by="rownumber")

school_rands <- rbind(school_majorcities_s1_rand,school_majorcities_s2_rand,school_majorcities_s3_rand,
                      school_regremote_s1_rand,school_regremote_s2_rand,school_regremote_s3_rand) %>%
  dplyr::select(schoolid,treatment) %>%
  rename(LD = treatment)

in2_time1 <- left_join(in2_time1_pre,school_rands,by=c("schoolid2" = "schoolid")) %>%
  rename(LD = LD.y) %>%
  dplyr::select(-.imp,-.id,-LD.x) %>%
  right_join(pace_preds,by=c("tname","LD")) %>%
  filter(!is.na(schoolid)) %>%
  #if old LD is the same as new LD then keep the existing fu_pace value (even if missing)
  #if old LD is different to new LD, then pred becomes the fu_pace value
  mutate(fu_total_pa_fin = case_when(is.na(fu_total_pa_fin) ~ NA,
                                     old_LD == LD ~ fu_total_pa_fin,
                                     old_LD != LD ~ pred))

in2_time1_split <- split(x=in2_time1,f=in2_time1$.imp)

#model the first interim
in2.time1 <- brms::brm_multiple(fu_total_pa_fin | mi() ~ LD + ses + RA_name + stype + base_total_pa_fin + trial + (1| schoolid), 
                                data = in2_time1_split,cores=4,chains=4,
                                iter = 5000, control = list(adapt_delta = 0.95), 
                                verbose = T,  save_pars = save_pars(all=TRUE), seed = 26085054,
                                prior = c(prior(normal(0,1000),class=b),
                                          prior(normal(0,1000),class=Intercept),
                                          prior(student_t(3, 0, 62), lb=0,class=sd),
                                          prior(student_t(3, 0, 62), lb=0,class=sigma)))

in2.time1.samp <- brms::posterior_samples(in2.time1)[,2:3] %>%
  data.frame() %>%
  #Want to see how many times HD compared to LD > 16.4 difference
  mutate(NI_margin = case_when(b_LDL - b_LDH >= -16.4 ~ 1, #LD is not inferior to HD by more than 16.4 mins
                               b_LDL - b_LDH < -16.4 ~ 0), #LD IS inferior to HD by more than 16.4 mins
         LD_HD_sup = case_when(b_LDL - b_LDH > 0 ~ 1, #LD is superior to HD
                               b_LDL - b_LDH <= 0 ~ 0), #LD is inferior to HD
         HD_sup_margin = case_when(b_LDH > 0 ~ 1,
                                   b_LDH <= 0 ~ 0),
         LD_sup_margin = case_when(b_LDL > 0 ~ 1,
                                   b_LDL <= 0 ~ 0))%>%
  summarise(NI_margin = mean(NI_margin),
            LD_HD_sup = mean(LD_HD_sup),
            HD_sup_margin = mean(HD_sup_margin),
            LD_sup_margin = mean(LD_sup_margin))

avg_comparisons(in2.time1,variables=list(LD="pairwise"),allow_new_levels=TRUE,cross=TRUE,ndraws=1000)
#Conclusion: Pr(noninferior) = 0.99
#early stopping for efficacy
#NI_margin
#LD_HD_sup
#HD_sup_margin
#LD_sup_margin
#0.991337
#0.836526
#0.999881
#0.999984
#Estimate 2.5 % 97.5 %                C: LD
#44.0  23.1   64.0 mean(H) - mean(Ctrl)
#55.5  30.8   87.0 mean(L) - mean(Ctrl)
#12.7 -11.9   40.3 mean(L) - mean(H)   
