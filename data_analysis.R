#### Data analysis ####
library(tidyverse)

# load data
load(file = "~/SEM/data/nhanes_cov_1516.RData")
load(file = "~/SEM/data/nhanes_out_1516.RData")
load(file = "~/SEM/data/nhanes_metals_1516.RData")
# load(file = "~/SEM/data/nhanes_phalates_pfas_1516.RData")

# log trasform chemicals
df_metals_log = df_metals %>% 
  select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_metals$SEQN, .) %>%
  mutate(SEQN = df_metals$SEQN) %>%
  select(- "df_metals$SEQN")

# join data
df_out_analysis = df_out %>% select(SEQN,BPXSY1,BPXDI1,
                                    BMXWAIST,BMXBMI)
df = join_all(list(df_cov,
                   df_out_analysis,
                   df_metals_log), 
              by='SEQN', type='full')
df = df %>% select(-SEQN)
dim(df)


# maybe we want to log transform some of the outcomes and some of the covariates
