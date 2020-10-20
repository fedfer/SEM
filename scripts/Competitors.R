# source competitors fuctions
library(tidyverse)
source("./scripts/Competitors_fcts_uni.R")

# read the data 
Y <- readRDS("~/SEM/results/Y.rds")
X_imputed <- readRDS("~/SEM/results/X_imputed.rds") %>% 
  select(-c("URDTNE6", "URDTNE7",)) %>%
  as.matrix()

# BMI
ind_ok <- is.finite(Y$BMXBMI)
bmi <- Y$BMXBMI[ind_ok] %>% scale() %>% as.vector() 
X_bmi <-  X_imputed[ind_ok,] %>% scale() %>% as.matrix()
hiernet_out <- Hiernet_fct(bmi, X_bmi)
family_out <- FAMILY_fct(bmi, X_bmi)
RAMP_out <- RAMP_fct(bmi, X_bmi)
list_bmi <- list(hiernet = hiernet_out, 
                 family = family_out, 
                 ramp = RAMP_out)
saveRDS(list_bmi, "./results/bmi_competitors.rds")
