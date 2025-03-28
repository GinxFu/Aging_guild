library(openxlsx)

# read metadata 
data_metadata <- read.xlsx('metadata.xlsx')
data_metadata_healthy <- read.xlsx('metadata_healthy.xlsx')
data_metadata_unhealthy <- read.xlsx('metadata_unhealthy.xlsx')

# read microbial age data
data_Micro_age <- read.xlsx('Microbial age.xlsx')
 
# merge
data_Micro_age_Dis <- merge(data_metadata, data_Micro_age[,-2], by = 'ID', all.x = T)

# -----Difference in microbial age-----
Dis_group <- c('Diabetes_h', 'Obesity_h', 'FattyDis_h', 'Gout_h', 
               'Hypertension_h', 'Stroke_h', 'HeartDis_h', 'Malnutrition_h')
results <- as.data.frame(matrix(NA, nrow = length(Dis_group), ncol = 9,
                                dimnames = list(c('T2DM', 'Obesity', 'FattyDis', 'Gout', 
                                                  'Hypertension', 'Stroke', 'HeartDis', 'Underweight'),
                                                c('Est','Se','t','P',
                                                  'Est_adj','Se_adj','t_adj','P_adj',
                                                  'Interaction_P'))))
for (i in 1:length(Dis_group)) {
  data_loop <- data_Micro_age_Dis[data_Micro_age_Dis$Sex == 0, ] # Sex == 0 for women or Sex == 1 for men
  names(data_loop)[names(data_loop) == Dis_group[i]] <- 'Dis'
  data_loop <- data_loop[data_loop$Dis %in% c('0','1'), ]
  
  model_lm.loop <- lm(Micro_age ~ Dis + Age, 
                      data = data_loop)
  
  results[i,c(1:4)] <- summary(model_lm.loop)$coefficients[2,]
  
  if (i %in% c(2,8)) {
    model_lm.loop <- lm(Micro_age ~ Dis + 
                          Age + Energy_kcal + 
                          smoking_7d + education + PA_frequently + Base_drk 
                        , data = data_loop)
  } else {
    model_lm.loop <- lm(Micro_age ~ Dis + 
                          Age + BMI + Energy_kcal + 
                          smoking_7d + education + PA_frequently + Base_drk 
                        , data = data_loop)
  }
  
  results[i,c(5:8)] <- summary(model_lm.loop)$coefficients[2,]
  
  model_lm.loop <- lm(Micro_age ~ Dis + Age + Age*Dis, 
                      data = data_loop)
  results[i,c(9)] <- summary(model_lm.loop)$coefficients[4,4]
  
  rm(data_loop, model_lm.loop)
}

# significance level
results <- within(results, {
  P_adj_star <- NA
  P_adj_star [P_adj < 0.10] <- '#'
  P_adj_star [P_adj < 0.05] <- '*'
  P_adj_star [P_adj < 0.01] <- '**'
  P_adj_star [P_adj < 0.001] <- '***'
  
  P_star <- NA
  P_star [P < 0.10] <- '#'
  P_star [P < 0.05] <- '*'
  P_star [P < 0.01] <- '**'
  P_star [P < 0.001] <- '***'
})

results_men <- results
results_women <- results

write.xlsx(results_men, 'MA_Dis_men.xlsx', rowNames = T)
write.xlsx(results_women, 'MA_Dis_women.xlsx', rowNames = T)



# -----T2DM risk-----
library(tidycmprsk)

# Outcome data
data_T2DM_cmp <- data_Micro_age_Dis[data_Micro_age_Dis$Rlt_healthy == '1', ] # only in healthy subjects
data_T2DM_cmp$Outcome_follow <- data_T2DM_cmp$DM_follow
data_T2DM_cmp$Outcome_follow [data_T2DM_cmp$Death_follow == 1] <- 2
table(data_T2DM_cmp$Outcome_follow)

# Old and Young microbial age
data_T2DM_cmp$Micro_age_diff <- data_T2DM_cmp$Micro_age -  data_T2DM_cmp$fit
data_T2DM_cmp <- within(data_T2DM_cmp, {
  Micro_age_diff_group <- NA
  Micro_age_diff_group [Micro_age_diff < 0] <- 'Young'
  Micro_age_diff_group [Micro_age_diff > 0] <- 'Old'
})
data_T2DM_cmp$Micro_age_diff_group <- factor(data_T2DM_cmp$Micro_age_diff_group, levels = c('Young', 'Old'))

# prepare data for analysis
data_analysis <- data_T2DM_cmp[data_T2DM_cmp$Sex == 1,] # men
data_analysis <- data_T2DM_cmp[data_T2DM_cmp$Sex == 0,] # women
data_analysis <- data_T2DM_cmp # all

table(data_analysis$Micro_age_diff_group, data_analysis$Outcome_follow)


# F-G test
cuminc.model <- cuminc(Surv(Months, as.factor(Outcome_follow)) ~ Micro_age_diff_group, data = data_analysis)
cuminc.model
sprintf("%0.3f", cuminc.model$cmprsk$Tests[1,2])


# CMP cox regression
cmp.model <- crr(Surv(Months, as.factor(Outcome_follow)) ~ Micro_age_diff_group 
                 # + Sex # only for all subjects
                 + Age + BMI + Energy_kcal
                 + education + PA + smoking + drinking, 
                 data = data_analysis)
cmp.model

HR <- sprintf("%0.2f", exp(c(cmp.model$tidy$estimate[1], cmp.model$tidy$conf.low[1], cmp.model$tidy$conf.high[1]))) 
HR
paste0(HR[1], " (", HR[2], ", ", HR[3], ")")
sprintf("%0.3f", cmp.model$tidy$p.value[1])


# CMP cox regression for difference 
cmp.model_diff <- crr(Surv(Months, as.factor(Outcome_follow)) ~ Micro_age_diff 
                      # + Sex # only for all subjects
                      + Age + BMI + Energy_kcal
                      + education + PA + smoking + drinking, 
                      data = data_analysis)
cmp.model_diff

HR_diff <- sprintf("%0.2f", exp(c(cmp.model_diff$tidy$estimate[1], cmp.model_diff$tidy$conf.low[1], cmp.model_diff$tidy$conf.high[1]))) 
HR_diff
paste0(HR_diff[1], " (", HR_diff[2], ", ", HR_diff[3], ")")
sprintf("%0.3f", cmp.model_diff$tidy$p.value[1])

