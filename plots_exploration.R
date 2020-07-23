# Plots from Infinite factor 
library(tidyverse)
library(plyr)
library(dplyr)
library(mice)
library(Rcpp)
library(RcppArmadillo)
library(infinitefactor)

# load data (local)--------------------------
load(file = "data/nhanes_cov_1516.RData")
load(file = "data/nhanes_out_1516.RData")
load(file = "data/nhanes_chem_1516.RData")

# # Log transform chemicals--------------------------
# TODO: maybe we want to log transform some of the outcomes and some of the covariates
df_chem_log = df_chem %>% 
  dplyr::select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_chem$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_chem$SEQN) %>% #so we append new column named SEQN to the end of dataset
  dplyr::select(- "df_chem$SEQN") #and delete column with ugly name

df_chem_log <- df_chem_log %>%
  dplyr::select(-c("LBXBCO", "LBXBCR", "LBXCOT", "LBXHCT", "URXANBT", "URXCOTT", "URXHCTT", "URDTNE2", "URXANTT", "URXNICT", "URXNNCT", "LBDWFL", "URX4MDA", "URX4TDA", "URX5NDA", "URX6TDA", "URXPPDA", "LBDSF1LC", "LBDSF2LC", "LBDSF3LC", "LBDSF4LC", "LBDSF5LC", "LBDSF6LC", "LBXIHG", "LBXBCD", "LBXBGE", "LBXBGM", "LBXTHG", "URXUCD", "URXUHG", "URXUMN", "URXUSB", "URXUTU", "LBXMPAH", "LBXPFDE", "LBXPFHS", "SSACET", "SSAND", "SSCLOT", "SSIMID", "SSOHIM", "SSTHIA", "URXUUR", "LBXBFOA", "LBXPFDO", "LBXPFUA", "URXMC1", "URXMCOH", "URXMHBP", "URXMHNC", "URXMHP", "LBXTST", "URXMNP", "URXUAB", "URXUAC", "URXUAS3", "URXUAS5", "URXUDMA", "URXUMMA", "LBX2DF", "LBX4CE", "LBXEST", "LBXV06", "LBXV07N", "LBXV08N", "LBXV3B", "LBXV4C", "URXUTRI"))

df_chem_log <- df_chem_log %>%
  dplyr::select(-(LBXVBF:LBXVDB))

df_chem_log <- df_chem_log %>%
  dplyr::select(-(LBXVDE:LBXVMC))

df_chem_log <- df_chem_log %>%
  dplyr::select(-(LBXVMCP:LBXVTP))

df_chem_log <- df_chem_log %>%
  dplyr::select(-c("LBXVVB", "LBXVXY", "URXBPM", "URXCYHA", "URXCYM", "URX1DC", "URX2DC", "URXDPM", "URXGAM", "URXIPM1", "URXIPM3", "URXMB1", "URXMB2", "URXPHE", "URXPMA", "URXTCV"))

# Examine covariates ------------------
df_cov <- df_cov %>% dplyr::select(-c("DMDEDUC3", "DMDEDUC2", "DMDMARTL", "RIDEXPRG", "LBXAPB", "URXUAS", "URXUCL")) # remove covariates with too many missing values


# Outcomes------------------------
df_out_analysis = df_out %>% dplyr::select(SEQN,BPXSY1,BPXDI1,
                                           BMXWAIST,BMXBMI)


# Join dataset--------------------------
df = join_all(list(df_out_analysis,
                   df_cov,
                   df_chem_log), 
              by='SEQN', type='full')
df = df %>% dplyr::select(-SEQN)

df_out_complete = df_out_analysis[complete.cases(df_out_analysis),] %>% 
  dplyr::select(-SEQN)

# Try linear Factor Model
out = linearMGSP(X = df_out_complete %>% scale() %>% as.matrix(), nrun = 1000, burn = 500,
                 adapt = FALSE, kinit = 2)
aligned = jointRot(out$lambdaSamps, out$etaSamps)
plotmat(lmean(aligned$lambda))

eig_val = eigen(cor(df_out_complete))$values
cumsum(eig_val)/sum(eig_val)




\