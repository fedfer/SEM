library(ggplot2)
library(infinitefactor)

# gibbs_out <- readRDS(file = "results/gibbs_results_Fede.rds")  # I didn't store loadings matrix and latent variables in this run
gibbs_out <- readRDS(file = "gibbs_results.rds")
str(gibbs_out)
plot(gibbs_out$acp) # Problem with storing acp in gibbs sampler

iterations <- dim(gibbs_out$Phi_st)[1]

chem_names <- read.table("chem_names_test.txt")
chem_names <- as.vector(chem_names[,1])


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


vec_true_inter_coeff_est_BMI <- c(true_inter_coeff_est_BMI) # vectorize

plot(vec_true_inter_coeff_est_BMI)


# 95% Credible intervals for interactions coefficients-------

# lower intervals
inter_coeff_lower_interval <- apply(gibbs_out$inter_coeff_st, c(2, 3, 4), function(x){
  interval <- quantile(x, probs = c(0.025, 0.975))
  return(interval[1])
})
inter_coeff_lower_interval_BMI <- inter_coeff_lower_interval[4,,]
# add up symmetric lower intervals
inter_coeff_lower_interval_BMI_upper_tri <- matrix(data = 0, nrow = nrow(inter_coeff_lower_interval_BMI), 
                                                   ncol = ncol(inter_coeff_lower_interval_BMI))
for (i in 1:nrow(inter_coeff_lower_interval_BMI_upper_tri)) {
  for (j in i:ncol(inter_coeff_lower_interval_BMI_upper_tri)) {
    if (i != j) {
      inter_coeff_lower_interval_BMI_upper_tri[i, j] <- inter_coeff_lower_interval_BMI[i, j] + inter_coeff_lower_interval_BMI[j, i]
    }
    else {
      inter_coeff_lower_interval_BMI_upper_tri[i, j] <- inter_coeff_lower_interval_BMI[i, j]
    }
  }
}

# upper intervals
inter_coeff_upper_interval <- apply(gibbs_out$inter_coeff_st, c(2, 3, 4), function(x){
  interval <- quantile(x, probs = c(0.025, 0.975))
  return(interval[2])
})
inter_coeff_upper_interval_BMI <- inter_coeff_upper_interval[4,,]
# add up symmetric upper intervals
inter_coeff_upper_interval_BMI_upper_tri <- matrix(data = 0, nrow = nrow(inter_coeff_upper_interval_BMI), 
                                                   ncol = ncol(inter_coeff_upper_interval_BMI))
for (i in 1:nrow(inter_coeff_upper_interval_BMI_upper_tri)) {
  for (j in i:ncol(inter_coeff_upper_interval_BMI_upper_tri)) {
    if (i != j) {
      inter_coeff_upper_interval_BMI_upper_tri[i, j] <- inter_coeff_upper_interval_BMI[i, j] + inter_coeff_upper_interval_BMI[j, i]
    }
    else {
      inter_coeff_upper_interval_BMI_upper_tri[i, j] <- inter_coeff_upper_interval_BMI[i, j]
    }
  }
}

# sanity check on the intervals
sum((inter_coeff_upper_interval_BMI - inter_coeff_lower_interval_BMI) < 0) # passed
sum((inter_coeff_upper_interval_BMI_upper_tri - inter_coeff_lower_interval_BMI_upper_tri) < 0) # passed


# interaction coefficients with 0 not in credible interval-----
tmp_nrow = nrow(true_inter_coeff_est_BMI)
tmp_ncol = ncol(true_inter_coeff_est_BMI)
mat_select_BMI <- matrix(data = 0, nrow = tmp_nrow, ncol = tmp_ncol) # Quick way to select entrie using zero one or true false matrix?
for (i in 1:tmp_nrow) {
  for (j in i:tmp_ncol) {
    if (0 < inter_coeff_lower_interval_BMI_upper_tri[i, j] | 0 > inter_coeff_upper_interval_BMI_upper_tri[i, j]) {
      mat_select_BMI[i, j] = 1
    }
  }
}
sum(mat_select_BMI) #number of interactions selected

signif_inter_BMI <- matrix(data = NA, nrow = sum(mat_select_BMI), ncol = 5)
colnames(signif_inter_BMI) <- c("chem_inter_name", "coeff", "lower_interval", "upper_interval", "method")

# examine upper triangular part
count <- 1

for (i in 1:tmp_nrow) {
  for (j in i:tmp_ncol) {
    if (mat_select_BMI[i, j] == 1) {
      signif_inter_BMI[count, "chem_inter_name"] <- paste(chem_names[i], "x", chem_names[j])
      signif_inter_BMI[count, "coeff"] <- true_inter_coeff_est_BMI[i, j]
      signif_inter_BMI[count, "lower_interval"] <- inter_coeff_lower_interval_BMI_upper_tri[i, j]
      signif_inter_BMI[count, "upper_interval"] <- inter_coeff_upper_interval_BMI_upper_tri[i, j]
      signif_inter_BMI[count, "method"] <- "SEM" 
      count <- count + 1
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
signif_main_BMI <- matrix(data = NA, nrow = sum(vec_select_BMI), ncol = 5)
colnames(signif_main_BMI) <- c("chem_name", "coeff", "lower_interval", "upper_interval", "method")
#rownames_signif_main_BMI <- vector(mode = "character", length = nrow(signif_main_BMI))
j <- 1
for (i in 1:tmp_length) {
  if (vec_select_BMI[i] == 1) {
    signif_main_BMI[j, "coeff"] <- true_coeff_est_BMI[i]
    signif_main_BMI[j, "chem_name"] <- chem_names[i]
    signif_main_BMI[j, "lower_interval"] <- main_coeff_lower_interval_BMI[i]
    signif_main_BMI[j, "upper_interval"] <- main_coeff_upper_interval_BMI[i]
    j <- j + 1
  }
}
signif_main_BMI[, "method"] <- rep("SEM", times = nrow(signif_main_BMI))



#convert to df-----
#convert main to df
df_signif_main_BMI <- as.data.frame(signif_main_BMI)
df_signif_main_BMI[, 2:4] <- lapply(df_signif_main_BMI[, 2:4], function(x) as.numeric(as.character(x)))

#convert interactions to df
df_signif_inter_BMI <- as.data.frame(signif_inter_BMI)
df_signif_inter_BMI[, 2:4] <- lapply(df_signif_inter_BMI[, 2:4], function(x) as.numeric(as.character(x)))



# Get chemical names for interactions-----------


# Visualizations-----
ggplot(df_signif_main_BMI, aes(x=chem_name, y=coeff)) + 
  geom_pointrange(aes(ymin=lower_interval, ymax=upper_interval)) +
  scale_y_continuous(breaks = round(seq(min(df_signif_main_BMI$lower_interval), max(df_signif_main_BMI$upper_interval), by = 0.005),2)) +
  theme(axis.text.x = element_text(angle=45))

ggplot(df_signif_inter_BMI, aes(x=chem_inter_name, y=coeff)) + 
  geom_pointrange(aes(ymin=lower_interval, ymax=upper_interval)) +
  scale_y_continuous(breaks = round(seq(min(df_signif_inter_BMI$lower_interval), max(df_signif_inter_BMI$upper_interval), by = 0.005),2)) +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle=90))
#fix the zero-one selection matrix for interaction elements


