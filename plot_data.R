# Load libraries and data
library(tidyverse)
library(plyr)
library(reshape2)
load("~/SEM/nhanes_1516.RData")


# matrix 0-1 for missing data
df_01 = df %>% is.na()
df_01 %>% as.matrix() %>% image()


# Correlation matrix
# Correlation plot
Cor_plot = cor(df %>% as.matrix(),use = "pairwise.complete.obs")
colnames(Cor_plot) = rownames(Cor_plot) = colnames(df)
Cor_plot = melt(Cor_plot)
ggplot(Cor_plot, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value), colour="grey20") + 
  scale_fill_gradient2(low = "#191970", high = "#006400", mid = "white") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_text(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(fill = " ") + 
  ggtitle("Correlation matrix")



# Plot two matrices
library(reshape2)
SampleMean = cbind(melt(sample_mean), label)
label = "Aligned Sample Mean"
ProcessMean = Reduce("+", aligned)/length(aligned)
rownames(ProcessMean) = colnames(X)
ProcessMean = cbind(melt(ProcessMean), label)
ggdf = rbind(SampleMean, ProcessMean)
ggplot(ggdf, aes(x = Var2, y = Var1)) + 
  facet_grid(cols = vars(label)) +
  geom_tile(aes(fill=value), colour="grey20") + 
  scale_fill_gradient2(low = "#800000", high = "#006400", mid = "white") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        #axis.text = element_blank(),
        legend.title = element_text(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(fill = " ")
