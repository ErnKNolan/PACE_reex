
pacman::p_load(mice,brms,here,tidyverse,gtsummary)
'%!in%' <- function(x,y)!('%in%'(x,y))

pace <- readRDS(here("Data","pace.RDS")) %>%
  mutate(pacetrial = factor(ifelse(trial == 0,"Noninferiority","Superiority"),
                            levels=c("Superiority","Noninferiority"))) %>%
  group_by(schoolid) %>%
  mutate(base_pa = mean(base_total_pa_fin,na.rm=T))

demographics_school <- pace %>% 
  group_by(schoolid) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>%
  tbl_strata(strata = pacetrial, ~.x %>%
             tbl_summary(by = LD,include = c("stype","ses","RA_name","base_pa")))

demographics_schooldf <- data.frame(demographics_school)

demo_diffs <- pace %>% group_by(schoolid) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  tbl_summary(by=pacetrial, include = c("stype","ses","RA_name","base_pa")) %>%
  add_p()

demo_diffsdf <- data.frame(demo_diffs)
