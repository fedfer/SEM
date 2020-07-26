library(infinitefactor)
library(coda)
library(lattice)

chem_names <- scan(file = "chem_names_test.txt", what = character()) #read in names of chemicals to use later

gibbs_out <- readRDS(file = "gibbs_results_2.rds") # Loading matrices included

iterations <- dim(gibbs_out$Phi_st)[1]



# Examine the loadings matrices lambda using infinitefactor package-----

# loading matrix for the outcome Y
Lambda_y_est <- apply(gibbs_out$Lambda_y_st, c(2, 3), mean) # first get posterior mean as estimate
plotmat(varimax(Lambda_y_est)[[1]])

# loading matrix for the covariates X
Lambda_x_est <- apply(gibbs_out$Lambda_x_st, c(2, 3), mean) # first get posterior mean as estimate
plotmat(varimax(Lambda_x_est)[[1]])
gibbs_out$
#resolve rotational ambituity for loading matrix of Y
# didn't store xi and eta!! gonna rerun



# Diagnostics-----
coeff_st_BMI <- gibbs_out$coeff_st[,4,]
colnames(coeff_st_BMI) <- chem_names
coeff_st_BMI_mcmc <- as.mcmc(coeff_st_BMI)
xyplot(coeff_st_BMI_mcmc[, 1:3]) # traceplot
densityplot(coeff_st_BMI_mcmc[, 1:3]) # check posterior distributions and assessing convergence
levelplot(coeff_st_BMI_mcmc) # correlation between the parameters
acfplot(coeff_st_BMI_mcmc[, 1:3]) #auto-correlation within chains, although we only have one chain
summary(coeff_st_BMI_mcmc) # summary statistics
effectiveSize(coeff_st_BMI_mcmc) # effective sample size




