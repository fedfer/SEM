#### Data analysis ####
library(tidyverse)
library(plyr)
library(dplyr)
library(mice)
library(Rcpp)
library(RcppArmadillo)

# load data (for server)--------------------------
load(file = "/work/sta790/ff31/SEM/data/nhanes_cov_1516.RData")
load(file = "/work/sta790/ff31/SEM/data/nhanes_out_1516.RData")
load(file = "/work/sta790/ff31/SEM/data/nhanes_chem_1516.RData")

# load data (local)--------------------------
# load(file = "data/nhanes_cov_1516.RData")
# load(file = "data/nhanes_out_1516.RData")
# load(file = "data/nhanes_chem_1516.RData")

# source (server)--------------------------
source("/work/sta790/ff31/SEM/scripts/gibbs_cpp_Fede.R")
 
# source (local)--------------------------
# source("scripts/gibbs_inter_cov_missing_v1.R")
# source("scripts/gibbs_cpp_Fede.R")

# # Log transform chemicals--------------------------
# TODO: maybe we want to log transform some of the outcomes and some of the covariates
df_chem_log = df_chem %>% 
  dplyr::select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_chem$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_chem$SEQN) %>% #so we append new column named SEQN to the end of dataset
  dplyr::select(- "df_chem$SEQN") #and delete column with ugly name

# Plot histogram (comment out for server)--------------------------
# df_chem_log[,1:10] %>%
#   gather() %>%
#   ggplot(aes(value)) +
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram()


# Throw away explanatory variables that are not normally distributed (metal)--------------------------
# variables thrown away: "LBXBCD", "LBXBCR", "LBXBGE", "LBXBGM", "LBXIHG", "LBXTHG", "URXUAB",
# "URXUAC", "URXUAS3", "URXUAS5", "URXUCD", "URXUDMA", "URXUHG", "URXUMMA", "URXUMN"
# df_metals = df_metals %>%
#   dplyr::select(- c("LBXBCD", "LBXBCR", "LBXBGE", "LBXBGM", "LBXIHG",
#                 "LBXTHG", "URXUAB", "URXUAC", "URXUAS3", "URXUAS5", "URXUCD", 
#                 "URXUDMA", "URXUHG", "URXUMMA", "URXUMN"))
# 
# df_metals_log = df_metals_log %>%
#   dplyr::select(- c("LBXBCD", "LBXBCR", "LBXBGE", "LBXBGM", "LBXIHG",
#               "LBXTHG", "URXUAB", "URXUAC", "URXUAS3", "URXUAS5", "URXUCD", 
#               "URXUDMA", "URXUHG", "URXUMMA", "URXUMN"))


# Throw away explanatory variables that are not normally distributed (chem)--------------------------
# Variables thrown away: "LBXBCO", "LBXBCR", "LBXCOT", "LBXHCT", "URXANBT", 
#                         "URXCOTT", "URXHCTT", "URDTNE2", "URXANTT", "URXNICT", 
#                         "URXNNCT", "LBDWFL", "URX4MDA", "URX4TDA", "URX5NDA", 
#                         "URX6TDA", "URXPPDA", "LBDSF1LC", "LBDSF2LC", "LBDSF3LC", 
#                       "LBDSF4LC", "LBDSF5LC", "LBDSF6LC", "LBXIHG", "LBXBCD", "LBXBGE", 
#                           "LBXBGM", "LBXTHG", "URXUCD", "URXUHG", "URXUMN", "URXUSB", "URXUTU",
#                             "LBXMPAH", "LBXPFDE", "LBXPFHS", "SSACET", "SSAND", "SSCLOT", 
#                           "SSIMID", "SSOHIM", "SSTHIA", "URXUUR", "LBXBFOA", "LBXPFDO", 
#                             "LBXPFUA", "URXMC1", "URXMCOH", "URXMHBP", "URXMHNC", "URXMHP", 
#                             "LBXTST", "URXMNP", "URXUAB", "URXUAC", "URXUAS3", "URXUAS5", "URXUDMA", "URXUMMA", 
#                             "LBX2DF", "LBX4CE", "LBXEST", "LBXV06", "LBXV07N", "LBXV08N",
#                               "LBXV3B", "LBXV4C", "URXUTRI"
# LBXVBF:LBXVDB
# LBXVDE:LBXVMC
# LBXVMCP:LBXVTP
# "LBXVVB", "LBXVXY", "URXBPM", "URXCYHA", "URXCYM", "URX1DC", "URX2DC", "URXDPM", "URXGAM", "URXIPM1", "URXIPM3", "URXMB1", "URXMB2", "URXPHE", "URXPMA", "URXTCV"


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
# summary(df_cov)
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

# Imputing covariates------------------------------------------
imp <- mice(df[, ncol(df_out_analysis):ncol(df)]) # when SEQN is included in df_out_analysis
imp <- complete(imp)


# Split into Y, X, Z--------------------------
Y <- df[, 1:4]
Z <- imp[, 1:(1+ncol(df_cov)-1-1)]
Z <- as.matrix(Z)
X <- df[,(ncol(Y) + ncol(Z) + 1):ncol(df)]

# Center Y and X ----
X_scaled <- scale(X, scale = TRUE)
Y_scaled <- scale(Y, scale = TRUE)
X <- X_scaled
Y <- Y_scaled

# 0-1 matrices for missing values, entry is 1 if missing--------------------------
Y_NA <- Y %>% is.na() * 1
X_NA <- X %>% is.na() * 1


# 0-1 matrices for lod, log scale (new)--------------------------------------------
chem_lod_log =  X %>% apply(., 2, function(x) min(x, na.rm = T))
chem_lod_log = chem_lod_log %>%
  rep(., times = nrow(X)) %>%
  matrix(., nrow = nrow(X), byrow = T)
chem_lod_matrix_log_scale = (chem_lod_log == X) %>% as.matrix() *1
chem_lod_matrix_log_scale[is.na(chem_lod_matrix_log_scale)] <- 0
X_LOD <- chem_lod_matrix_log_scale
LOD_X_vec <- chem_lod_log

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


X_hollow <- X_hollow - chem_lod_matrix_log_scale * (as.data.frame(X_hollow) %>% as.matrix())


# Sanity check: 0's in Y_hollow and X_hollow should conform to the 0-1 matrices ---------------------
colSums(X_hollow == 0) - colSums( (X_LOD + X_NA) == 1)
colSums(X == 0 & !is.na(X)) # Sanity check for X_hollow passed

# how many latent factors?------------------------
df_chem_imp <- imp[, (1+ncol(df_cov)-1):ncol(imp)]
dim(df_chem_imp)
cor_chem = df_chem_imp %>%
  cor(. ,use="complete.obs")
eig_chem = cor_chem %>% eigen()
plot(eig_chem$values)
vec <- cumsum(eig_chem$values)/sum(eig_chem$values)
sum(vec > 0.9)

cor(df_out_analysis, use="complete.obs")
cor_out = df_out_analysis %>% 
  dplyr::select(-SEQN) %>% 
  cor(. ,use="complete.obs") 
eig_out = cor_out %>% eigen()
plot(eig_out$values)
vec <- cumsum(eig_out$values)/sum(eig_out$values)
sum(vec > 0.9)



# Gibbs ------------------------------------------
nrun = 1000
burn = 500
n_samples = nrun - burn
gibbs_result <- gibbs(X = X_hollow, Y = Y_hollow,
                      X_NA = X_NA, Y_NA = Y_NA, X_LOD = X_LOD, LOD_X_vec = LOD_X_vec, Z = Z,
                      nrun = nrun, burn = burn, thin = 1, alpha_prior = NULL, theta_inf = 0.05,
                      k = 32, m = 2, a = 1/2, delta_rw = 0.5)

results_dir="/work/sta790/ff31/SEM/"
saveRDS(gibbs_result, file.path(results_dir, "gibbs_results_Fede.rds"))









