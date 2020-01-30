#### Data analysis ####
library(tidyverse)
library(plyr)
library(dplyr)
library(mice)

# load data--------------------------
load(file = "/work/yj90/SEM/data/nhanes_cov_1516.RData")
load(file = "/work/yj90/SEM/data/nhanes_out_1516.RData")
load(file = "/work/yj90/SEM/data/nhanes_metals_1516.RData")
load(file = "/work/yj90/SEM/data/nhanes_cov_simple_1516.RData")

# source--------------------------
source("/work/yj90/SEM/scripts/gibbs_inter_cov_missing_v1.R")

# # Log transform chemicals--------------------------
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


# Examine covariates ------------------
summary(df_cov)
df_cov <- df_cov %>% select(-c("DMDEDUC3", "DMDEDUC2", "DMDMARTL", "RIDEXPRG", "LBXAPB", "URXUAS", "URXUCL")) # remove covariates with too many missing values


# Outcomes------------------------
df_out_analysis = df_out %>% dplyr::select(SEQN,BPXSY1,BPXDI1,
                                           BMXWAIST,BMXBMI)

# Plot histogram outcomes--------------------------
df_out_analysis[,] %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
# BMXBMI looks a bit skewed, otherwise looks okay


# Join dataset--------------------------

df = join_all(list(df_out_analysis,
                   df_cov,
                   df_metals_log), 
              by='SEQN', type='full')

df = df %>% dplyr::select(-SEQN)

# Imputing covariates------------------------------------------
imp <- mice(df[, ncol(df_out_analysis):ncol(df)]) # when SEQN is included in df_out_analysis
imp <- complete(imp)


# Split into Y, X, Z--------------------------
Y <- df[, 1:4]
Z <- imp[, 1:(1+ncol(df_cov)-1-1)]
X <- df[,(ncol(Y) + ncol(Z) + 1):ncol(df)]


# 0-1 matrices for missing values, entry is 1 if missing--------------------------
Y_NA <- Y %>% is.na() * 1
X_NA <- X %>% is.na() * 1

# 0-1 matrices for lod, log scale --------------------------------------------
# metal_log_scale <- df %>% dplyr::select(LBXBCO:URXUUR) %>% as.matrix()
# metal_lod_log =  metal_log_scale %>% apply(., 2, function(x) min(x, na.rm = T))
# metal_lod_log = metal_lod_log %>%
#   rep(., times = nrow(metal_log_scale)) %>%
#   matrix(., nrow = nrow(metal_log_scale), byrow = T)
# metal_lod_matrix_log_scale = (metal_lod_log == metal_log_scale) %>% as.matrix() *1
# metal_lod_matrix_log_scale[is.na(metal_lod_matrix_log_scale)] <- 0
# colSums(metal_lod_matrix_log_scale)
# 
# X_LOD <- cbind( matrix(data = 0, nrow = nrow(df %>% dplyr::select(RIDAGEYR:LBXTC) %>% as.matrix()), ncol = ncol(df %>% dplyr::select(RIDAGEYR:LBXTC) %>% as.matrix())),
#                 metal_lod_matrix_log_scale )
# 
# 
# LOD_X_vec <- c( rep(0, times = ncol(df %>% dplyr::select(RIDAGEYR:LBXTC) %>% as.matrix())),
#                 metal_log_scale %>% apply(., 2, function(x) min(x, na.rm = T)))

# 0-1 matrices for lod, log scale (new)--------------------------------------------
metal_lod_log =  X %>% apply(., 2, function(x) min(x, na.rm = T))
metal_lod_log = metal_lod_log %>%
  rep(., times = nrow(X)) %>%
  matrix(., nrow = nrow(X), byrow = T)
metal_lod_matrix_log_scale = (metal_lod_log == X) %>% as.matrix() *1
metal_lod_matrix_log_scale[is.na(metal_lod_matrix_log_scale)] <- 0
X_LOD <- metal_lod_matrix_log_scale
LOD_X_vec <- metal_lod_log

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
# X_hollow <- X_hollow - ( cbind( matrix(data = 0, nrow = nrow(df %>% dplyr::select(RIDAGEYR:LBXTC) %>% as.matrix()), ncol = ncol(df %>% dplyr::select(RIDAGEYR:LBXTC) %>% as.matrix())),
#                                  metal_lod_matrix_log_scale * (as.data.frame(X_hollow) %>% dplyr::select(LBXBCO:LBXTC) %>% as.matrix()) ) )

X_hollow <- X_hollow - metal_lod_matrix_log_scale * (as.data.frame(X_hollow) %>% as.matrix())


# Sanity check: 0's in Y_hollow and X_hollow should conform to the 0-1 matrices ---------------------
colSums(X_hollow == 0) - colSums( (X_LOD + X_NA) == 1)
colSums(df == 0 & !is.na(df)) # Sanity check for X_hollow passed

# how many latent factors?------------------------
dim(df_metals_log)
cor_metals = df_metals_log %>% 
  dplyr::select(-SEQN) %>% 
  cor(. ,use="complete.obs") 
eig_metals = cor_metals %>% eigen()
plot(eig_metals$values)
cumsum(eig_metals$values)/sum(eig_metals$values)


cor(df_out_analysis, use="complete.obs")
cor_out = df_out_analysis %>% 
  dplyr::select(-SEQN) %>% 
  cor(. ,use="complete.obs") 
eig_out = cor_out %>% eigen()
plot(eig_out$values)
cumsum(eig_out$values)/sum(eig_out$values)



# Gibbs ------------------------------------------
nrun = 1000
burn = 500
n_samples = nrun - burn
gibbs_result <- gibbs(X = X_hollow, Y = Y_hollow,
                      X_NA = X_NA, Y_NA = Y_NA, X_LOD = X_LOD, LOD_X_vec = LOD_X_vec, Z = Z_test,
                      nrun = nrun, burn = burn, thin = 1, alpha_prior = NULL, theta_inf = 0.05,
                      k = 7, m = 3, a = 1/2, delta_rw = 0.1) # maybe startoff delta_rw as 0.05, acceptance rate high at first then became too low for 0.1









