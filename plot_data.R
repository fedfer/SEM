# Load libraries and data
library(tidyverse)
library(plyr)
library(reshape2)
library(naniar)


load("~/SEM/data/nhanes_1516.RData")
load(file = "data/nhanes_metals_1516.RData")

# Plot histogram Metals
df_metals_log = df_metals %>% 
  select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_metals$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_metals$SEQN) %>% #so we append new column named SEQN to the end of dataset
  select(- c("df_metals$SEQN","SEQN"))

df_metals_log[,11:22] %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

quantile(df_metals$LBXBCR, na.rm = T)
quantile(df_metals$LBXIHG, na.rm = T)

# matrix 0-1 for missing data
load("data/nhanes_1516.RData")
df_01 = df %>% is.na()
df_01 %>% as.matrix() %>% image()
#vis_miss(df_chem, warn_large_data = F)


# matrix 0-1 for Limit of Detection
#vis_miss(df_out, warn_large_data = F)
load("data/nhanes_chem_1516.RData")
df_chem_noSEQN = df_chem %>% dplyr::select(-SEQN)
lod =  df_chem_noSEQN %>% apply(., 2, function(x) min(x, na.rm = T))
chem_lod = df_chem_noSEQN %>% apply(., 2, function(x) min(x, na.rm = T)) %>%
  rep(., times = nrow(df_chem_noSEQN)) %>%
  matrix(., nrow = nrow(df_chem_noSEQN), byrow = T)

df_lod = (df_chem_noSEQN == chem_lod) %>% as.matrix() 
df_lod %>% as.matrix() %>% image()
lod_sum = apply(df_lod, 2 , function(x)  sum(x,na.rm = T))
hist(lod_sum)

# Correlation plot outcomes
load("data/nhanes_out_1516.RData")
Cor_plot = cor(df_out %>% as.matrix(),use = "pairwise.complete.obs")
colnames(Cor_plot) = rownames(Cor_plot) = colnames(df_out)
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
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_text(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(fill = " ") + 
  ggtitle("Correlation Outcomes")

# Correlation plot chemicals
Cor_plot = cor(df_chem %>% as.matrix(),use = "pairwise.complete.obs")
colnames(Cor_plot) = rownames(Cor_plot) = colnames(df_chem)
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
  ggtitle("Correlation Chemicals")



# Plot two matrices
library(reshape2)
df_01_melt = melt(df_01)
ggplot(df_01_melt, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value), colour="grey20") + 
  #scale_fill_gradient2(low = "#800000", high = "#006400", mid = "white") +
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
