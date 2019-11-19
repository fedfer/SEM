#### Data analysis ####
library(tidyverse)
library(plyr)
library(dplyr)

# load data--------------------------
load(file = "data/nhanes_cov_1516.RData")
load(file = "data/nhanes_out_1516.RData")
load(file = "data/nhanes_metals_1516.RData")

# source ----------------------------------------
source("scripts/gibbs_interactions_missing_v1.R")

# Log transform chemicals--------------------------
# TODO: maybe we want to log transform some of the outcomes and some of the covariates
df_metals_log = df_metals %>% 
  dplyr::select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_metals$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_metals$SEQN) %>% #so we append new column named SEQN to the end of dataset
  dplyr::select(- "df_metals$SEQN") #and delete column with ugly name

# Plot histogram metals--------------------------
df_metals_log[,] %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()


# Throw away explanatory variables that are not normally distributed--------------------------
# variables thrown away: "LBXBCD", "LBXBCR", "LBXBGE", "LBXBGM", "LBXIHG", "LBXTHG", "URXUAB",
# "URXUAC", "URXUAS3", "URXUAS5", "URXUCD", "URXUDMA", "URXUHG", "URXUMMA", "URXUMN"
df_metals = df_metals %>%
  dplyr::select(- c("LBXBCD", "LBXBCR", "LBXBGE", "LBXBGM", "LBXIHG", "LBXTHG", "URXUAB", "URXUAC", "URXUAS3", "URXUAS5", "URXUCD", "URXUDMA", "URXUHG", "URXUMMA", "URXUMN"))

df_metals_log = df_metals_log %>%
  dplyr::select(- c("LBXBCD", "LBXBCR", "LBXBGE", "LBXBGM", "LBXIHG", "LBXTHG", "URXUAB", "URXUAC", "URXUAS3", "URXUAS5", "URXUCD", "URXUDMA", "URXUHG", "URXUMMA", "URXUMN"))


# Get values limit of detection in metals (original scale)--------------------------
df_metals_noSEQN <- df_metals %>% dplyr::select(-SEQN)
lod <- df_metals_noSEQN %>% apply(., 2, function(x){min(x, na.rm = TRUE)})

# Join dataset--------------------------
df_out_analysis = df_out %>% dplyr::select(SEQN,BPXSY1,BPXDI1,
                                    BMXWAIST,BMXBMI)
df = join_all(list(df_out_analysis,
                   df_cov,
                   df_metals_log), 
              by='SEQN', type='full')

df = df %>% dplyr::select(-SEQN)

# Plot histogram outcomes--------------------------
df_out_analysis[,] %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
# BMXBMI looks a bit skewed, otherwise looks okay

# modify lod to suite dimension of joined dataset--------------------------
# lod <- c(rep(0, times = ncol(df_cov) - 1 ), lod) # length of this vector should be the same as the number of columns of X that will come later


# Split into Y and X--------------------------
Y <- df[, 1:4]
X <- df %>% dplyr::select(-c("BPXSY1", "BPXDI1", "BMXWAIST", "BMXBMI"))

# 0-1 matrices for missing values, entry is 1 if missing--------------------------
Y_NA <- Y %>% is.na() * 1
X_NA <- X %>% is.na() * 1

# 0-1 matrices for lod --------------------------------------------
df_orig_scale <- join_all(list(df_out_analysis,
                               df_cov,
                               df_metals), 
                          by='SEQN', type='full')
df_orig_scale <- df_orig_scale %>% dplyr::select(-SEQN)

metal_orig_scale <- df_orig_scale %>% dplyr::select(LBXBCO:URXUUR) %>% as.matrix()

metal_lod <- lod %>% rep(., times = nrow(metal_orig_scale)) %>%
  matrix(., nrow = nrow(metal_orig_scale), byrow = T)

metal_lod_matrix_orig_scale <- (metal_lod == metal_orig_scale) %>% as.matrix() *1
metal_lod_matrix_orig_scale[is.na(metal_lod_matrix_orig_scale)] <- 0
colSums(metal_lod_matrix_orig_scale)

# Compare results with Federico's code, in log scale
metal_log_scale <- df %>% dplyr::select(LBXBCO:URXUUR) %>% as.matrix()
metal_lod_log =  metal_log_scale %>% apply(., 2, function(x) min(x, na.rm = T))
metal_lod_log = metal_lod_log %>%
  rep(., times = nrow(metal_log_scale)) %>%
  matrix(., nrow = nrow(metal_log_scale), byrow = T)
metal_lod_matrix_log_scale = (metal_lod_log == metal_log_scale) %>% as.matrix() *1
metal_lod_matrix_log_scale[is.na(metal_lod_matrix_log_scale)] <- 0
colSums(metal_lod_matrix_log_scale)

X_LOD <- cbind( matrix(data = 0, nrow = nrow(df %>% dplyr::select(RIDAGEYR:URXUCL) %>% as.matrix()), ncol = ncol(df %>% dplyr::select(RIDAGEYR:URXUCL) %>% as.matrix())),
                metal_lod_matrix_log_scale )


LOD_X_vec <- c( rep(0, times = ncol(df %>% dplyr::select(RIDAGEYR:URXUCL) %>% as.matrix())),
                metal_log_scale %>% apply(., 2, function(x) min(x, na.rm = T)))



# Set places where there is NA to 0 in X and Y-----------------------
Y_hollow <- apply(Y, c(1, 2), function(x){
  if (is.na(x)) {
    return(0)
  }else{
    return(x)
  }
})

X_hollow <- apply(X, c(1, 2), function(x){
  if (is.na(x)) {
    return(0)
  }else{
    return(x)
  }
})

# Set places where there is LOD to 0 ---------------------
# metals
X_hollow <- X_hollow - ( cbind( matrix(data = 0, nrow = nrow(df %>% dplyr::select(RIDAGEYR:URXUCL) %>% as.matrix()), ncol = ncol(df %>% dplyr::select(RIDAGEYR:URXUCL) %>% as.matrix())),
                                metal_lod_matrix_log_scale * (as.data.frame(X_hollow) %>% dplyr::select(LBXBCO:URXUUR) %>% as.matrix()) ) )


# Sanity check: 0's in Y_hollow and X_hollow should conform to the 0-1 matrices ---------------------
colSums(X_hollow == 0) - colSums( (X_LOD + X_NA) == 1)
colSums(df == 0 & !is.na(df)) # Sanity check for X_hollow passed

# Gibbs ------------------------------------------
nrun = 1000
burn = 500
n_samples = nrun - burn
gibbs_result <- gibbs(X = X_hollow, Y = Y_hollow, X_NA = X_NA, Y_NA = Y_NA, X_LOD = X_LOD, LOD_X_vec = LOD_X_vec,
                      nrun = nrun, burn = burn, thin = 1,
                      alpha_prior = NULL, theta_inf = 0.05,
                      m = ncol(Y_hollow), k = ncol(X_hollow))

# Debug------------------------------------------
# Got NA in my Psi!
# sum(is.na(X_hollow))
# sum(is.na(Y_hollow))
# sum(is.na(X_NA))
# sum(is.na(Y_NA))
# sum(is.na(X_LOD))
# sum(is.na(LOD_X_vec))



 



