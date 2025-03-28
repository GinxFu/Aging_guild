
# read metadata 
data_metadata_healthy <- read.xlsx('metadata_healthy.xlsx')
data_metadata_unhealthy <- read.xlsx('metadata_unhealthy.xlsx')

# read microbial age data
data_Micro_age_men <- read.xlsx('Microbial age men.xlsx')
data_Micro_age_women <- read.xlsx('Microbial age women.xlsx')
data_Micro_age <- rbind(data_Micro_age_men, data_Micro_age_women); data_Micro_age <- data_Micro_age[order(data_Micro_age$ID), ]; rownames(data_Micro_age) <- NULL
rm(data_Micro_age_men, data_Micro_age_women)

# merge with metadata
data_Micro_age_healthy <- data_Micro_age[data_Micro_age$ID %in% data_metadata_healthy$ID, ]
data_Micro_age_healthy <- merge(data_metadata_healthy, data_Micro_age_healthy[,-2], by = 'ID', all.x = T)
rm(data_Micro_age)


# Old and Young microbial age
data_Micro_age_healthy$Micro_age_diff <- data_Micro_age_healthy$Micro_age -  data_Micro_age_healthy$fit
data_Micro_age_healthy <- within(data_Micro_age_healthy, {
  Micro_age_diff_group <- NA
  Micro_age_diff_group [Micro_age_diff < 0] <- 'Young'
  Micro_age_diff_group [Micro_age_diff > 0] <- 'Old'
})
data_Micro_age_healthy$Micro_age_diff_group <- factor(data_Micro_age_healthy$Micro_age_diff_group, levels = c('Young', 'Old'))

table(data_Micro_age_healthy$Micro_age_diff_group, useNA = 'ifany')
table(data_Micro_age_healthy$Sex, data_Micro_age_healthy$Micro_age_diff_group, useNA = 'ifany')


# LE analysis
library(flexsurv)
library(rms)
library(survival)
library(pracma)

model.flexsurv <- flexsurvspline(Surv(Months, Death_follow) ~ rcs(Age, 4)  + Sex + Micro_age_diff_group, 
                                 k = 4, data = data_Micro_age_healthy)


# compute LE
seq <- seq(45,100,0.1)
result_LE <- data.frame(age = seq)
result_LE$surv_prob_1O <- result_LE$surv_prob_1Y <- result_LE$surv_prob_0O <- result_LE$surv_prob_0Y <- NA
for (i in 1:nrow(result_LE)) {
  
  # women Young
  pred_prob <- predict(model.flexsurv, type = 'surv',newdata = data.frame(Age = result_LE$age[i], Sex = 0,
                                                                          Micro_age_diff_group = 'Young'))[[1]][[1]]
  surv_prob <- mean(pred_prob$.pred_survival)
  result_LE$surv_prob_0Y[i] <- surv_prob
  
  # women Old
  pred_prob <- predict(model.flexsurv, type = 'surv',newdata = data.frame(Age = result_LE$age[i], Sex = 0,
                                                                          Micro_age_diff_group = 'Old'))[[1]][[1]]
  surv_prob <- mean(pred_prob$.pred_survival)
  result_LE$surv_prob_0O[i] <- surv_prob
  
  # men Young
  pred_prob <- predict(model.flexsurv, type = 'surv',newdata = data.frame(Age = result_LE$age[i], Sex = 1,
                                                                          Micro_age_diff_group = 'Young'))[[1]][[1]]
  surv_prob <- mean(pred_prob$.pred_survival)
  result_LE$surv_prob_1Y[i] <- surv_prob
  
  # men Old
  pred_prob <- predict(model.flexsurv, type = 'surv',newdata = data.frame(Age = result_LE$age[i], Sex = 1,
                                                                          Micro_age_diff_group = 'Old'))[[1]][[1]]
  surv_prob <- mean(pred_prob$.pred_survival)
  result_LE$surv_prob_1O[i] <- surv_prob

  rm(pred_prob, surv_prob)
}


# integrating
data_loop <- result_LE
result_LE$LE_1O <- result_LE$LE_1Y <- result_LE$LE_0O <- result_LE$LE_0Y <- NA

for (i in 1:nrow(result_LE)) {
  if (i > 1) {
    data_loop <- data_loop[-1,]
  }
  result_LE$LE_0Y[i] <- trapz(data_loop$age, data_loop$surv_prob_0Y)
  result_LE$LE_0O[i] <- trapz(data_loop$age, data_loop$surv_prob_0O)
  result_LE$LE_1Y[i] <- trapz(data_loop$age, data_loop$surv_prob_1Y)
  result_LE$LE_1O[i] <- trapz(data_loop$age, data_loop$surv_prob_1O)
}

write.xlsx(result_LE, 'Life expectance.xlsx')
