# source competitors fuctions
library(tidyverse)
source("./scripts/Competitors_fcts_uni.R")

# read the data 
Y <- readRDS("~/SEM/results/Y.rds")
X_imputed <- readRDS("~/SEM/results/X_imputed.rds") %>% 
  select(-c("URDTNE6", "URDTNE7")) %>%
  as.matrix()

# BMI
bmi <- Y$BMXBMI %>% as.vector()
ind_ok <- is.finite(bmi)
hiernet_out <- Hiernet_fct(bmi[ind_ok],X_imputed[ind_ok,])
family_out <- FAMILY_fct(bmi[ind_ok],X_imputed[ind_ok,])
RAMP_out <- RAMP_fct(bmi[ind_ok],X_imputed[ind_ok,])
list_bmi <- list(hiernet = hiernet_out, 
                 family = family_out, 
                 ramp = RAMP_out)
saveRDS(list_bmi, "./results/bmi_competitors.rds")


# systolic 
syst <- Y$BPXSY1 %>% as.vector()
ind_ok <- is.finite(syst)
hiernet_out <- Hiernet_fct(syst[ind_ok],X_imputed[ind_ok,])
family_out <- FAMILY_fct(syst[ind_ok],X_imputed[ind_ok,])
RAMP_out <- RAMP_fct(syst[ind_ok],X_imputed[ind_ok,])
list_syst <- list(hiernet_out, family_out, family_out)
saveRDS(list_syst, "./results/syst_competitors.rds")