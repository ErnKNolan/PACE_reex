imp <- readRDS(here("Data","imp.RDS"))

pace_hds <- imp %>% filter(LD == "H") %>%
  mutate(trial = factor(ifelse(as.numeric(as.character(schoolid)) %in% c(
    2,4,6,10,11,14,15,16,19,21,23,25,26,27,31,32,34,37,38,39,42,43,47,
    49,50,51,52,54,55,60,61
  ),"Sup","Noninf")))

vs_hds <- brms::brm(fu_total_pa_fin | mi() ~ trial + ses + RA_name + stype + base_total_pa_fin + (1| schoolid), 
                    data = pace_hds, 
                    iter = 2500, control = list(adapt_delta = 0.90), 
                    verbose = T,  save_pars = save_pars(all=TRUE), seed = 26085054,
                    cores=4)
summary(vs_hds)


vs_hds.samp <- brms::posterior_samples(vs_hds)[,2] %>%
  data.frame() %>%
  mutate(noninf_better = case_when(`.` < 0 ~ 1, 
                                   `.` >= 0 ~ 0))%>%
  summarise(noninf_better = mean(noninf_better))

avg_comparisons(vs_hds,variables=list(trial="pairwise"),allow_new_levels=TRUE,cross=TRUE,ndraws=1000)
