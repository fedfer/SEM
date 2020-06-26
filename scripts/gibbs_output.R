library(ggplot2)

gibbs_out <- readRDS(file = "gibbs_results.rds")

plot(gibbs_out$acp) # Problem with updating acp in gibbs sampler

iterations <- dim(gibbs_out$Phi_st)[1]

# par(mfrow=c(2,2))
# plot(x = 1:500, y = gibbs_out$Phi_st[, 1, 1])
# plot(x = 1:500, y = gibbs_out$Phi_st[, 2, 2])
# plot(x = 1:500, y = gibbs_out$Phi_st[, 3, 3])
# plot(x = 1:500, y = gibbs_out$Phi_st[, 4, 4])

# index <- 1
# for (i in 1:dim(gibbs_out$Phi_st)[2]) {
#   for (j in 1:dim(gibbs_out$Phi_st)[3]) {
#     p <- ggplot(x = 1:500, y = gibbs_out$Phi_st[, i, j])
#     assign(paste0('p', index), p)
#     index <- index + 1
#   }
# }
# rm(index)


par(mfrow=c(2,2))
for (i in 1:4) {
  plot(x = 1:500, y = gibbs_out$Phi_st[, i, i])
}

dev.off()

# compute coefficients for interactions
inter_coeff_est_st <- apply(gibbs_out$inter_coeff_st, c(2, 3, 4), mean) # What type of value is this?
# analyze coefficients for BMI
inter_coeff_est_BMI <- inter_coeff_est_st[4,,]

# compute actual coefficients
p <- dim(inter_coeff_est_BMI)[1]
true_inter_coeff_est_BMI <- matrix(data = 0, nrow = p, ncol = p)

for (i in 1:p) {
  for (j in 1:p) {
    if(i != j){
      true_inter_coeff_est_BMI[i, j] <- inter_coeff_est_BMI[i, j] + inter_coeff_est_BMI[j, i]
    }
    else {
      true_inter_coeff_est_BMI[i, j] <- inter_coeff_est_BMI[i, j]
    }
  }
}

heatmap(true_inter_coeff_est_BMI)

vec_true_inter_coeff_est_BMI <- c(true_inter_coeff_est_BMI)

par(mfrow=c(2,1))

plot(vec_true_inter_coeff_est_BMI)

# values of coefficients greater than 0.05
inter_coeff_est_BMI_0.05 <- true_inter_coeff_est_BMI
inter_coeff_est_BMI_0.05[abs(inter_coeff_est_BMI_0.05) < 0.05]=0 
plot(c(inter_coeff_est_BMI_0.05))

dev.off()

# values of coefficients greater than 0.1                  # how small counts?
inter_coeff_est_BMI_0.1 <- true_inter_coeff_est_BMI
inter_coeff_est_BMI_0.1[abs(inter_coeff_est_BMI_0.1) < 0.1]=0 
plot(c(inter_coeff_est_BMI_0.1))

# Credible intervals-------


# Main effects-------------------
true_coeff_est <- apply(gibbs_out$coeff_st, c(2, 3), mean)
true_coeff_est_BMI <- true_coeff_est[4,]
plot(true_coeff_est_BMI) # coefficients for main effects seems to be smaller than interactions


# Get chemical names for main effects------------
coeff_est_BMI_0.05 <- true_coeff_est_BMI
coeff_est_BMI_0.05[abs(true_coeff_est_BMI) < 0.05]=0
plot(coeff_est_BMI_0.05)
chem_names[coeff_est_BMI_0.05 != 0]

# Get chemical names for interactions-----------

