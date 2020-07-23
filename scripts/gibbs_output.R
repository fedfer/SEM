library(ggplot2)
library(infinitefactor)

gibbs_out <- readRDS(file = "gibbs_results_Feder.rds")  # I didn't store loadings matrix and latent variables in this run

plot(gibbs_out$acp) # Problem with storing acp in gibbs sampler

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
inter_coeff_est_st <- apply(gibbs_out$inter_coeff_st, c(2, 3, 4), mean)
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

# heatmap(true_inter_coeff_est_BMI)
# dev.off()

vec_true_inter_coeff_est_BMI <- c(true_inter_coeff_est_BMI) # vectorize


plot(vec_true_inter_coeff_est_BMI)

# values of coefficients greater than 0.05
# inter_coeff_est_BMI_0.05 <- true_inter_coeff_est_BMI
# inter_coeff_est_BMI_0.05[abs(inter_coeff_est_BMI_0.05) < 0.05]=0 
# plot(c(inter_coeff_est_BMI_0.05))
# 
# dev.off()

# values of coefficients greater than 0.1                 
# inter_coeff_est_BMI_0.1 <- true_inter_coeff_est_BMI
# inter_coeff_est_BMI_0.1[abs(inter_coeff_est_BMI_0.1) < 0.1]=0 
# plot(c(inter_coeff_est_BMI_0.1))

# 95% Credible intervals for interactions coefficients-------
inter_coeff_lower_interval <- apply(gibbs_out$inter_coeff_st, c(2, 3, 4), function(x){
  interval <- quantile(x, probs = c(0.025, 0.975))
  return(interval[1])
})
inter_coeff_lower_interval_BMI <- inter_coeff_lower_interval[4,,]

inter_coeff_upper_interval <- apply(gibbs_out$inter_coeff_st, c(2, 3, 4), function(x){
  interval <- quantile(x, probs = c(0.025, 0.975))
  return(interval[2])
})
inter_coeff_upper_interval_BMI <- inter_coeff_upper_interval[4,,]

# sanity check on the intervals
sum((inter_coeff_upper_interval_BMI - inter_coeff_lower_interval_BMI) < 0) # passed

# interaction coefficients with 0 not in credible interval-----
tmp_nrow = nrow(true_inter_coeff_est_BMI)
tmp_ncol = ncol(true_inter_coeff_est_BMI)
mat_select_BMI <- matrix(data = 0, nrow = tmp_nrow, ncol = tmp_ncol) # Quick way to select entrie using zero one or true false matrix?
for (i in 1:tmp_nrow) {
  for (j in 1:tmp_ncol) {
    if (0 < inter_coeff_lower_interval_BMI[i, j] | 0 > inter_coeff_upper_interval_BMI[i, j]) {
      mat_select_BMI[i, j] = 1
    }
  }
}
signif_inter_BMI <- matrix(data = NA, nrow = sum(mat_select_BMI), ncol = 1)
for (i in 1:tmp_nrow) {
  for (i in 1:tmp_ncol) {
    if (mat_select_BMI[i, j] == 1) {
      ### How to deal with symmetry here when considering credible intervals for interactions? ###
    }
  }
}


# Main effects-------------------
true_coeff_est <- apply(gibbs_out$coeff_st, c(2, 3), mean)
true_coeff_est_BMI <- true_coeff_est[4,]

plot(c(true_coeff_est_BMI))

# 95% Credible intervals for main effects-------
main_coeff_lower_interval <- apply(gibbs_out$coeff_st, c(2, 3), function(x){
  interval <- quantile(x, probs = c(0.025, 0.975))
  return(interval[1])
})
main_coeff_lower_interval_BMI <- main_coeff_lower_interval[4,]

main_coeff_upper_interval <- apply(gibbs_out$coeff_st, c(2, 3), function(x){
  interval <- quantile(x, probs = c(0.025, 0.975))
  return(interval[2])
})
main_coeff_upper_interval_BMI <- main_coeff_upper_interval[4,]

# main coefficients with 0 not in credible interval-----
tmp_length <- length(true_coeff_est_BMI)
vec_select_BMI <- rep(0, times = tmp_length)
for (i in 1:tmp_length) {
  if (0 < main_coeff_lower_interval_BMI[i] | 0 > main_coeff_upper_interval_BMI[i] ) {
    vec_select_BMI[i] = 1
  }
}
signif_main_BMI <- matrix(data = NA, nrow = sum(vec_select_BMI), ncol = 1)
rownames_signif_main_BMI <- vector(mode = "character", length = nrow(signif_main_BMI))
j <- 1
for (i in 1:tmp_length) {
  if (vec_select_BMI[i] == 1) {
    signif_main_BMI[j, 1] = true_coeff_est_BMI[i]
    rownames_signif_main_BMI[j] <- chem_names[i]
    j <- j + 1
  }
}
rownames(signif_main_BMI) <- rownames_signif_main_BMI
signif_main_BMI




# Get chemical names for main effects------------
coeff_est_BMI_0.05 <- true_coeff_est_BMI
coeff_est_BMI_0.05[abs(true_coeff_est_BMI) < 0.05]=0
plot(coeff_est_BMI_0.05)
chem_names[coeff_est_BMI_0.05 != 0]

# Get chemical names for interactions-----------

# Analysis for Lambda-----


