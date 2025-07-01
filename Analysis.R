#AUTHOR: Erin Nolan
#DATE CREATED: 18 Feb 2025
#PURPOSE: Analyse the data

pacman::p_load(mice,brms,here,tidyverse,gtsummary,marginaleffects)
'%!in%' <- function(x,y)!('%in%'(x,y))

pace <- readRDS(here("Data","pace.RDS"))

set.seed(866380771)
#set.seed(254782510) #sensitivity seed
pace_imp <- mice(pace,m=1,maxit=100,meth = c("","","","","","","pmm","",""),
                        pred=quickpred(pace,exclude=c("tname","fu_total_pa_fin")))
events <- pace_imp$loggedEvents
imp <- complete(pace_imp,action="long")
saveRDS(imp,here("Data","imp.RDS"))
#saveRDS(imp,here("Data","imp_sens.RDS"))
#split <- split(x=imp,f=imp$.imp)

#This is the full model
#impmodel <- bf(fu_total_pa_fin | mi() ~ LD*stype + ses + RA_name + base_total_pa_fin + trial + (1| schoolid))
impmodel <- bf(fu_total_pa_fin | mi() ~ LD + ses + RA_name + stype + base_total_pa_fin + trial + (1| schoolid))
#looks like 1 imputation for baseline pa is fine, doesnt change much and very little effect on the other vars
#10000 iterations for the main model but thats not needed for the later models
mod.adj <- brms::brm(impmodel, data = imp, 
                     chains=4,cores=4,
                     iter = 5000, control = list(adapt_delta = 0.95), 
                     verbose = T,  save_pars = save_pars(all=TRUE), #primary seed = 26085054, #sensitivity seed = 254782510
                     seed = 26085054,
                     prior = c(prior(normal(0,1000),class=b),
                               prior(normal(0,1000),class=Intercept),
                               prior(student_t(3, 0, 62), lb=0,class=sd),
                               prior(student_t(3, 0, 62), lb=0,class=sigma)))
summary(mod.adj)
#results
main.samp <- brms::posterior_samples(mod.adj)[,2:3] %>%
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
#NI_margin
#LD_HD_sup
#HD_sup_margin
#LD_sup_margin
#0.9802
#0.4739
#1
#1
avg_comparisons(mod.adj,variables=list(LD="pairwise"),allow_new_levels=TRUE,cross=TRUE,ndraws=10000)
#Estimate 2.5 % 97.5 %                C: LD
#47.283  36.2   58.6 mean(H) - mean(Ctrl)
#46.886  28.2   65.7 mean(L) - mean(Ctrl)
#-0.476 -15.7   14.4 mean(L) - mean(H)   
#avg_comparisons(mod.adj,variables=list(LD="pairwise",stype="all"),allow_new_levels=TRUE,cross=TRUE,ndraws=10000)#int effects
#rhats <- mod.adj$rhats

#the 3 counterfactuals
pace_LD <- imp %>% mutate(old_LD = LD,LD = "L")
pace_HD <- imp %>% mutate(old_LD = LD,LD = "H")
pace_ctrl <- imp %>% mutate(old_LD = LD,LD = "Ctrl")
pace_new <- rbind(pace_HD,pace_LD,pace_ctrl) %>%
  dplyr::select(-.imp)

#10 draws and 20, 30, 50 draws were very similar for treatment effects - keep at 25
preds <- brms::posterior_predict(mod.adj,newdata=pace_new,ndraws=100)

pace_list <- list()
for(i in 1:100){
  pace_list[[i]] <- data.frame(pace_new,pred=preds[i,],.imp=i)
}
pace_list2 <- bind_rows(pace_list)
pace_list2 <- pace_list2 %>%
  mutate(fu_total_pa_fin = case_when(old_LD == LD ~ fu_total_pa_fin,
                                     old_LD != LD ~ pred))
pace_preds <- pace_list2 %>% 
  dplyr::select(tname,LD,old_LD,pred,.imp)

saveRDS(pace_preds,here("Data","pace_preds_sens.RDS"))

#-------------------------------------------------------------
pace %>% group_by(schoolid) %>% filter(row_number()==1) %>%
tbl_summary(,by=trial,include=c(stype,LD,ses,RA_name,base_total_pa_fin)) %>%
  add_p()

#-------------------------------------------------------------

#Sensitivity analysis - not including the 2nd pace trial HD
pace <- readRDS(here("Data","pace.RDS"))
pace_sens <- pace %>% filter(as.numeric(as.character(schoolid)) %!in% c(
  2,4,6,10,11,14,15,16,19,21,23,25,26,27,31,32,34,37,38,39,42,43,47,
  49,50,51,52,54,55,60,61
)) %>%
  mutate(schoolid = factor(schoolid))
set.seed(866380771)
pace_sens_imp <- mice(pace_sens,m=1,maxit=100,meth = c("","","","","","","pmm","",""),
                 pred=quickpred(pace_sens,exclude=c("tname","fu_total_pa_fin")))
events <- pace_sens_imp$loggedEvents
sens_imp <- complete(pace_sens_imp,action="long")
saveRDS(sens_imp,here("Data","sens_imp.RDS"))

mod.adj.sens <- brms::brm(impmodel, data = sens_imp, 
                          iter = 10000, control = list(adapt_delta = 0.95), 
                          verbose = T,  save_pars = save_pars(all=TRUE), seed = 26085054,
                          prior = c(prior(normal(0,1000),class=b),
                                    prior(normal(0,1000),class=Intercept),
                                    prior(student_t(3, 0, 62), lb=0,class=sd),
                                    prior(student_t(3, 0, 62), lb=0,class=sigma)))
summary(mod.adj.sens)
#results
main.sens.samp <- brms::posterior_samples(mod.adj.sens)[,2:3] %>%
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

#the 3 counterfactuals
pace_LD <- sens_imp %>% mutate(old_LD = LD,LD = "L")
pace_HD <- sens_imp %>% mutate(old_LD = LD,LD = "H")
pace_ctrl <- sens_imp %>% mutate(old_LD = LD,LD = "Ctrl")
pace_new <- rbind(pace_HD,pace_LD,pace_ctrl) %>%
  dplyr::select(-.imp)
preds <- brms::posterior_predict(mod.adj.sens,newdata=pace_new,ndraws=10)

pace_list <- list()
for(i in 1:10){
  pace_list[[i]] <- data.frame(pace_new,pred=preds[i,],.imp=i)
}
pace_list2 <- bind_rows(pace_list)
pace_list2 <- pace_list2 %>%
  mutate(fu_total_pa_fin = case_when(old_LD == LD ~ fu_total_pa_fin,
                                     old_LD != LD ~ pred))
sens_pace_preds <- pace_list2 %>% 
  dplyr::select(tname,LD,old_LD,pred,.imp)

saveRDS(sens_pace_preds,here("Data","sens_pace_preds.RDS"))


#IF we do a simulation then 10000 draws, otherwise just do 1 draw of the counterfactual
#10000 draws
test <- brms::posterior_predict(mod.adj2,newdata=pace_new,ndraws=10000)
#posterior predict or posterior expected  - expected has same mean but smaller variance
#then based on which rep you're in, test[j,]
pace_new$pred <- test[2,]

#IF want to do a simulation part as well
#10000 reps for each combination
#but the rep is for the whole trial
#i.e you pull 1 rep, interim 1 you pull a draw, then interim 2 you pull another draw
#but you dont do 10000^10000 draws for 2 interims
#because the first block before interim 1 will always be the same, that stays the same
#but you then do multiple draws once you do counterfactuals


#schools in each arm / trial
#we would have to assume that instead of the 2 back to back, we instead had all arms
#start at the same time. So for HD start with sup HD then go into inf HD as trial progresses
#sup (HD) 31 schools, 2,4,6,10,11,14,15,16,19,21,23,25,26,27,31,32,34,37,38,39,42,43,47,49,50,51,52,54,55,60,61
#sup (ctrl) 30 schools 1,3,5,7,8,9,12,13,17,18,20,22,24,28,29,33,35,36,40,41,44,45,46,48,53,56,57,58,59,62
#inf (HD) 24 schools, 64,67,68,71,72,74,76,78,80,81,84,85,89,90,110,113,117,118,121,122,127,128,132,133
#inf (LD) 24 schools, 63,65,66,69,70,73,75,77,79,82,83,86,87,88,107,108,111,112,114,115,119,120,130,131

#can do some sensitivity seeing what would happen if HD was ordered differently, presented in a supplementary or something


#ICC getting
pace_s <- imp %>% filter(trial==1)
impmodel <- bf(fu_total_pa_fin | mi() ~ LD + ses + RA_name + stype + base_total_pa_fin + trial + (1| schoolid))
#looks like 1 imputation for baseline pa is fine, doesnt change much and very little effect on the other vars
#10000 iterations for the main model but thats not needed for the later models
mod.adj.s <- brms::brm(impmodel, data = pace_s, 
                     chains=4,cores=4,
                     iter = 5000, control = list(adapt_delta = 0.95), 
                     verbose = T,  save_pars = save_pars(all=TRUE), seed = 26085054,
                     prior = c(prior(normal(0,1000),class=b),
                               prior(normal(0,1000),class=Intercept),
                               prior(student_t(3, 0, 62), lb=0,class=sd),
                               prior(student_t(3, 0, 62), lb=0,class=sigma)))
summary(mod.adj.s)
performance::variance_decomposition(mod.adj.s)

pace_ni <- imp %>% filter(trial==0)
impmodel <- bf(fu_total_pa_fin | mi() ~ LD + ses + RA_name + stype + base_total_pa_fin + trial + (1| schoolid))
#looks like 1 imputation for baseline pa is fine, doesnt change much and very little effect on the other vars
#10000 iterations for the main model but thats not needed for the later models
mod.adj.ni <- brms::brm(impmodel, data = pace_ni, 
                       chains=4,cores=4,
                       iter = 5000, control = list(adapt_delta = 0.95), 
                       verbose = T,  save_pars = save_pars(all=TRUE), seed = 26085054,
                       prior = c(prior(normal(0,1000),class=b),
                                 prior(normal(0,1000),class=Intercept),
                                 prior(student_t(3, 0, 62), lb=0,class=sd),
                                 prior(student_t(3, 0, 62), lb=0,class=sigma)))
summary(mod.adj.ni)
performance::variance_decomposition(mod.adj.ni)