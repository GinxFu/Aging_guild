library(openxlsx)

# read microbial age data
data_Micro_age <- read.xlsx('Microbial age.xlsx')

# read metadata 
data_metadata_healthy <- read.xlsx('metadata_healthy.xlsx')

ID_men <- data_metadata_healthy$ID [data_metadata_healthy$Sex == 1]
ID_women <- data_metadata_healthy$ID [data_metadata_healthy$Sex == 0]


# -----define healthy lifestyles-----
data_analysis <- data_Pheno_healthy

## healthy diet
data_analysis <- within(data_analysis, {
  Healthy_diet <- NA
  Base_healthy_diet <- NA
  
  Healthy_diet [(Fiber > 20) | (CHFP > quantile(data_Pheno$CHFP, 0.60))] <- '1'
  Healthy_diet [! ((Fiber > 20) | (CHFP > quantile(data_Pheno$CHFP, 0.60)))] <- '0'
  
  Base_healthy_diet [(Fiber > 20) | (Base_CHFP > quantile(data_Pheno$Base_CHFP, 0.60))] <- '1'
  Base_healthy_diet [! ((Fiber > 20) | (Base_CHFP > quantile(data_Pheno$Base_CHFP, 0.60)))] <- '0'
})
      
## compute HLI
data_analysis <- within(data_analysis, {
  HLI_base <- NA
  HLI_follow <- NA
  HLI_long <- NA
  
  HLI_base <- 
    (Base_PA == 1) + (Base_smk == 0) + (Base_drk == '0') +
    (Base_healthy_diet == '1')
  
  HLI_follow <- 
    (PA == 1) + (smoking == 0) + (drinking == '0') +
    (Healthy_diet == '1')
  
  HLI_long [HLI_base <= 2 & HLI_follow <= 2] <- 'Low-Low'
  HLI_long [HLI_base <= 2 & HLI_follow > 2] <- 'Low-High'
  HLI_long [HLI_base > 2 & HLI_follow <= 2] <- 'High-Low'
  HLI_long [HLI_base > 2 & HLI_follow > 2] <- 'High-High'
  
})
data_analysis$HLI_long <- factor(data_analysis$HLI_long,
                                 levels = c('Low-Low',
                                            'High-Low',
                                            'Low-High',
                                            'High-High'))

# -----Analysis-----
data_analysis2 <- cbind(data_Micro_age[, -c(1,2)], data_analysis)
data_analysis2$Micro_age_diff <- data_analysis2$Micro_age - data_analysis2$fit


# linear regression

# Behaviors: Healthy_diet PA smoking drinking
# HLI: HLI_follow HLI_long

# data: data_analysis2
# data: data_analysis2[data_analysis2$ID %in% ID_men, ]
# data: data_analysis2[data_analysis2$ID %in% ID_women, ]


model.lm <- lm(Micro_age ~ HL_score + Age + BMI + education 
               # + Sex # only for all subjects 
               , data = data_analysis2)
summary(model.lm)
