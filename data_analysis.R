#### Data analysis ####
library(tidyverse)
library(plyr)

# load data
load(file = "data/nhanes_cov_1516.RData")
load(file = "data/nhanes_out_1516.RData")
load(file = "data/nhanes_metals_1516.RData")
# load(file = "/data/nhanes_phalates_pfas_1516.RData")

# log trasform chemicals
df_metals_log = df_metals %>% 
  select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_metals$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_metals$SEQN) %>% #so we append new column named SEQN to the end of dataset
  select(- "df_metals$SEQN") #and delete column with ugly name

# join data
df_out_analysis = df_out %>% select(SEQN,BPXSY1,BPXDI1,
                                    BMXWAIST,BMXBMI)
df = join_all(list(df_cov,
                   df_out_analysis,
                   df_metals_log), 
              by='SEQN', type='full')
df = df %>% select(-SEQN)
dim(df) #sanity check: ncol(df_cov) + ncol(df_out_analysis) + ncol(df_metals_log) - 3 should match column numbers of df


# maybe we want to log transform some of the outcomes and some of the covariates






# Some exploratory stuffs --------------------------------------------------
summary(df)






# Experimenting with how different functions work------------------
play = df_metals %>% 
  select(-SEQN) %>% 
  log(., base = 10) %>% 
  mutate(SEQN = df_metals$SEQN)

  
  
  

  




