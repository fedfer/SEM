# Plots main effects and interactions

library(ggplot2)
library(infinitefactor)

gibbs <- readRDS(file = "results/gibbs_results.rds")
str(gibbs_out)

# BMI plots
bmi_comp <- readRDS(file = "results/bmi_competitors.rds")

# Plot main effects
names = read.table("chem_names_test.txt")
label = "MQR-SEM"
beta_hat = apply(gibbs$coeff_st[,4,], 2, mean)
q_sup = apply(gibbs$coeff_st[,4,], 2, function(x) quantile(x,probs = 0.95))
q_inf = apply(gibbs$coeff_st[,4,], 2, function(x) quantile(x,probs = 0.05))
gibbs_plot = data.frame(Values = beta_hat, 
                        Variables = names,
                        model = label,
                        q_inf = q_inf, q_sup = q_sup)
p <- length(beta_hat)
# label = "PIE"
# PIE_plot = data.frame(Values = PIE$beta[1:p], Variables = names,model = label, q_inf = PIE$beta[1:p],
#                       q_sup = PIE$beta[1:p])
label = "RAMP"
RAMP <- bmi_comp$ramp
RAMP_plot = data.frame(Values = RAMP$beta[1:p], Variables = names,
                       model = label, q_inf = RAMP$beta[1:p],
                       q_sup = RAMP$beta[1:p])
label = "Family"
Family <- bmi_comp$family
Family_plot = data.frame(Values = Family$beta[1:p], Variables = names,
                         model = label,q_inf = Family$beta[1:p],
                         q_sup = Family$beta[1:p])
label = "HierNet"
hiernet <- bmi_comp$hiernet
hiernet_plot = data.frame(Values = hiernet$beta[1:p], Variables = names,
                          model = label,q_inf = hiernet$beta[1:p],
                          q_sup = hiernet$beta[1:p])


# errors
ramp_err <- RAMP$err_pred ^.2 %>% mean()
family_err <- Family$err_pred ^.2 %>% mean()
hiernet_err <- hiernet$err_pred ^.2 %>% mean()
# y_hat = gibbs$y_pred %>% apply(. , 2, mean)
# (y - y_hat) %>% .^2 %>% mean()

# PLOT of MAin Effects
beta_plot = rbind(gibbs_plot,#PIE_plot,
                  #RAMP_plot,
                  Family_plot,hiernet_plot)
beta_plot <- beta_plot %>% mutate(Variables = V1)
ggplot(beta_plot, aes(x = Variables, y = Values, color = model,
                      shape = model))+
  geom_point(size = 2.3)+
  theme(axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        #panel.grid.major = element_blank(),
        panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        text = element_text(size=17),
        axis.text.x = element_text(angle = 90, hjust = 1),
        #legend.title = element_text(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(fill = " ") + 
  ylab("Estimated Main Effects") +
  # ggtitle("Estimated Main Effects")+
  geom_errorbar(aes(ymin=q_inf, ymax=q_sup), width=.1)


