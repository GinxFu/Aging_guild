
# read metadata 
data_metadata_healthy <- read.xlsx('metadata_healthy.xlsx')

ID_men <- data_metadata_healthy$ID [data_metadata_healthy$Sex == 1]
ID_women <- data_metadata_healthy$ID [data_metadata_healthy$Sex == 0]

# -----1. Alpha diversity-----
library(Maaslin2)

data_alpha <- read.xlsx('Alpha.xlsx')
data_alpha_healthy <- data_alpha[data_alpha$ID %in% data_metadata_healthy$ID, ]

data_analysis <- merge(data_metadata_healthy, data_alpha_healthy, by = 'ID')

# linear regression

# index: chao1, pielou
# data: data_analysis[ID_men, ], data_analysis[ID_women, ]
model.lm <- lm(index ~ Age + BMI + Energy_kcal
               # + education + PA + smoking + drinking # sensitive analysis
               , data = data_analysis)

results <- summary(model.lm)$coefficients

# interaction of sex*age

model.lm <- lm(index ~ Age + Sex + Age*Sex + BMI + Energy_kcal
               # + education + PA + smoking + drinking # sensitive analysis
               , data = data_analysis)

results_interaction <- summary(model.lm)$coefficients

# -----2. Beta diversity-----
library(vegan)

# 导入数据
data_beta <- read.csv('Bray-Curtis.csv', row.names = 1)
data_beta_healthy <- data_beta[data_metadata_healthy$ID, data_metadata_healthy$ID]


# Adonis/PERMANOVA
dist_matrix <- as.matrix(data_beta_healthy)
set.seed(123)
# dist: dist_matrix[ID_men, ID_men],  dist_matrix[ID_women, ID_women]
# data: data_metadata_healthy[ID_men, ],  data_metadata_healthy[ID_women, ]
RD <- adonis2(formula = dist_matrix ~ BMI + Energy_kcal 
              # + education + PA + smoking + drinking # sensitive analysis
              + Age, 
              data = data_metadata_healthy,
              permutaions = 999, method = "bray")
RD

# Adonis/PERMANOVA for interaction of sex*age
RD_interaction <- adonis2(formula = dist_matrix ~ BMI + Energy_kcal 
                          # + education + PA + smoking + drinking # sensitive analysis
                          + Age + Sex + Age*Sex, 
                          data = data_metadata_healthy,
                          permutaions = 999, method = "bray")
RD_interaction


# -----3. PCoA axis-----

# get PCOA information
PCOA <- cmdscale(dist_matrix, eig = T, k = 25) 
pcoa_point <- as.data.frame(PCOA$points) 
# merge with metadata
pcoa_data <- cbind(pcoa_point, data_metadata_healthy)  


results_PCoA25 <- as.data.frame(matrix(NA, nrow = 25, ncol = 10,
                                       dimnames = list(c(),
                                                       c('Axis','Var_exp','Var_exp_culm',
                                                         'Est_men', 'SE_men', 'P_men',
                                                         'Est_women', 'SE_women','P_women',
                                                         'P_interaction'))))
for (i in 1:25) {
  
  results_PCoA25$Axis [i] <- i
  # variation  explained
  results_PCoA25$Var_exp [i] <- round(PCOA2$eig[i] / sum(PCOA2$eig), 4) * 100  
  
  # cumulative variation  explained
  results_PCoA25$Var_exp_culm [i] <- sum(results_PCoA25$Var_exp [1:i])
  
  # linear regression with age
  data_loop <- pcoa2_data
  names(data_loop)[i] <- 'Axis'
  
  results_PCoA25$Est_men [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal, 
               data = data_loop[data_loop$Sex == 1, ]))$coefficients[2,1]
  results_PCoA25$SE_men [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal, 
               data = data_loop[data_loop$Sex == 1, ]))$coefficients[2,2]
  results_PCoA25$P_men [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal, 
               data = data_loop[data_loop$Sex == 1, ]))$coefficients[2,4]
  
  results_PCoA25$Est_women [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal, 
               data = data_loop[data_loop$Sex == 0, ]))$coefficients[2,1]
  results_PCoA25$SE_women [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal, 
               data = data_loop[data_loop$Sex == 0, ]))$coefficients[2,2]
  results_PCoA25$P_women [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal, 
               data = data_loop[data_loop$Sex == 0, ]))$coefficients[2,4]
  
  results_PCoA25$P_interaction [i] <- 
    summary(lm(Axis ~ Age + BMI + Energy_kcal + Age*Sex, 
               data = data_loop))$coefficients[6,4]
  
}

write.xlsx(results_PCoA25,'Bray_PCoA25.xlsx')


# -----4. Uniqueness-----

# compute uniqueness
data_uniqueness <- data.frame(Uniqueness = rep(NA,ncol(data_beta_healthy)))
rownames(data_uniqueness) <- rownames(data_beta_healthy)
for (i in 1:ncol(data_beta_healthy)) {
  data_uniqueness$Uniqueness[i] <- min(data_beta_healthy[-i,i])
}
data_uniqueness <- cbind(data_metadata_healthy, data_uniqueness)

# linear regression
# data: data_uniqueness[ID_men, ],  data_uniqueness[ID_women, ]
model_lm <- lm(Uniqueness ~ Age + BMI + Energy_kcal
               # + education + PA + smoking + drinking # sensitive analysis
               , data = data_uniqueness)
summary(model_lm)

# interaction of sex*age
model_lm_interaction <- lm(Uniqueness ~ Age + BMI + Energy_kcal + Age*Sex
                           # + education + PA + smoking + drinking # sensitive analysis
                           , data = data_uniqueness)
summary(model_lm_interaction)

# -----5. abundance-----
library(Maaslin2)

# read data of abundance
data_abundance <- read.csv('abundance_table.csv', row.names = 1)
data_abundance_healthy <- data_abundance[data_metadata_healthy$ID, ]

# clr transformation
for (i in 1:nrow(data_abundance_healthy)) {
  data_abundance_healthy[i,] [data_abundance_healthy[i,] == 0] <- 0.5
}
rm(i)
data_abundance_healthy <- data.frame(compositions::clr(data_abundance_healthy))

# MaAslin2
Model <- Maaslin2(
  input_data = data_abundance_healthy,  # data_abundance_healthy[ID_men, ],  data_abundance_healthy[ID_wpmen, ]
  input_metadata = data_metadata_healthy, # data_metadata_healthy[ID_men, ],  data_metadata_healthy[ID_women, ],
  output = 'Temp',
  min_abundance = -1000,
  min_prevalence = 0, 
  normalization = "NONE", 
  transform = "NONE", 
  fixed_effects = c('Age', 'BMI', 'Energy_kcal'),
  standardize = FALSE, 
  plot_heatmap = FALSE,
  plot_scatter = FALSE)

results <- Model$results
results <- results[results$metadata == 'Age', ]

# P value correction
results$pval_BH <- p.adjust(results$pval, method = 'BH')
results <- results[order(results$pval),]

write.xlsx(results, 'Abundance_men.xlsx')
write.xlsx(results, 'Abundance_women.xlsx')

# interaction of sex*age
# abundance_i is the colname of data_abundance
model_lm_interaction <- lm(abundance_i ~ Age + Sex + Age*Sex + BMI + Energy_kcal, 
                           data = cbind(data_metadata_healthy, data_abundance_healthy))
summary(model_lm_interaction)


