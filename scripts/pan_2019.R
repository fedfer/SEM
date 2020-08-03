library(Ball)
library(mice)
library(dplyr)

X <- readRDS(file = "X.rds")
X_imputed <- mice(X)
X_imputed <- complete(X_imputed)
saveRDS(X_imputed, file = "X_imputed.rds")

Y <- readRDS(file = "Y.rds")
Y_imputed <- mice(Y)
Y_imputed <- complete(Y_imputed)
saveRDS(Y_imputed, file = "Y_imputed.rds")

#mice refuses to impute URDTNE6 and URDTNE7
X_imputed <- X_imputed %>% dplyr::select(-c("URDTNE6", "URDTNE7"))


res <- bcorsis(y = Y_imputed, x = X_imputed) # doesn't handle missing data!

