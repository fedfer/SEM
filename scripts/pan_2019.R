library(Ball)
library(mice)
library(dplyr)

set.seed(42)

# X <- readRDS(file = "X.rds")
# X_imputed <- mice(X)
# X_imputed <- complete(X_imputed)
# saveRDS(X_imputed, file = "X_imputed.rds")
X_imputed <- readRDS(file = "/work/yj90/SEM/X_imputed.rds") # for server
# X_imputed <- readRDS(file = "X_imputed.rds")

# Y <- readRDS(file = "Y.rds")
# Y_imputed <- mice(Y)
# Y_imputed <- complete(Y_imputed)
# saveRDS(Y_imputed, file = "Y_imputed.rds")
Y_imputed <- readRDS(file = "/work/yj90/SEM/Y_imputed.rds") # for server
# Y_imputed <- readRDS(file = "Y_imputed.rds")

#mice refuses to impute URDTNE6 and URDTNE7
X_imputed <- X_imputed %>% dplyr::select(-c("URDTNE6", "URDTNE7"))


bcorsis_res <- bcorsis(y = Y_imputed, x = X_imputed) # doesn't handle missing data!

#locally main chemicals ranked by importance
bcorsis_res_chem_names <- chem_names[bcorsis_res[[1]]]
saveRDS(bcorsis_res_chem_names, file = "bcorsis_res_chem_names.rds")

results_dir="/work/yj90/SEM/" # for server
saveRDS(bcorsis_res, file.path(results_dir, "bcorsis_res.rds"))

