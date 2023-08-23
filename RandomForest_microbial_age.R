library(openxlsx)
library(randomForest)

# import data
Pheno <- read.xlsx('Pheno.xlsx', rowNames = T)
Pheno$Sex <- as.factor(Pheno$Sex)
Pheno$Rlt_health <- as.factor(Pheno$Rlt_health)
Agedata <- Pheno[c('Age', 'Sex', 'Rlt_health')]
CAG_table <- read.xlsx('CAG_table.xlsx', rowNames = T)

# clr-transformed
for (i in 1:130) {
  CAG_table[,i] [CAG_table[,i] == 0] <- 0.5
}
rm(i)
CAG_table <- data.frame(compositions::clr(CAG_table))

# keep healthy subjects and merge data
Agedata <- Agedata[Agedata$Rlt_health == '1',]
inmeta <- rownames(CAG_table) %in% rownames(Agedata)
CAGdata <- CAG_table[inmeta,]
rm(inmeta)
Age_CAG <- cbind(Agedata, CAGdata)   

# choose sex and analyze men and women separately
sex = '1'
# or
sex = '0'

Age_CAG <- Age_CAG[Age_CAG$Sex == sex, ]
Agedata <- Agedata[Agedata$Sex == sex, ]
CAGdata <- CAGdata[rownames(CAGdata) %in% rownames(Agedata), ]

# include significant CAG only
if (sex == '0') {
CAG_sig <- read.xlsx('Age_CAG_regression_women.xlsx')
CAG_sig <- CAG_sig[CAG_sig$pval < 0.05, 1]
} else if (sex == '1') {
CAG_sig <- read.xlsx('Age_CAG_regression_men.xlsx')
CAG_sig <- CAG_sig[CAG_sig$pval < 0.05, 1]
}

Age_CAG_train <- Age_CAG[c('Age', CAG_sig)]

# -----determine parameters-----
# cross-validation
set.seed(119)
rfcv.mod.re5 <- replicate(5, rfcv(CAGdata[CAG_sig], Agedata$Age, cv.fold = 10, step = 1.1), 
                          simplify = FALSE)
error.cv <- data.frame(sapply(rfcv.mod.re5, '[[', 'error.cv'))
error.cv$err.mean <- rowMeans(error.cv)
# error plot
matplot(rfcv.mod.re5[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type = "l",
        lwd = c(2, rep(1, ncol(error.cv))), col = 1, lty = 1, log = "x",
        xlab = "Number of CAGs", ylab = "CV Error")
# CV error reached a lowest level with a minimum of 20 features. 
# Therefore, we selected all 25 male-specific CAGs and 30 female-specific CAGs

# determine mtry
mtry.err <- data.frame(mtry = 0, err = 0)
for (i in 1:(ncol(Age_CAG_train) - 1)) {
  set.seed(119)
  mtry_fit <- randomForest(Age~., data = Age_CAG_train, mtry = i)
  err <- mean(mtry_fit$mse)
  mtry.err <- rbind(mtry.err, data.frame(mtry = i, err = err))
  rm(i, mtry_fit, err)
}
mtry.err <- mtry.err[-1,]; rownames(mtry.err) <- NULL
plot(mtry.err)
# include all 25 male-specific CAGs and 30 female-specific CAGs 
# mtry = 25 for men; mtry = 30 for women

# determine ntree
set.seed(119)
ntree_fit <- randomForest(Age~., data = Age_CAG_train, mtry = (ncol(Age_CAG_train) - 1), ntree = 2000)
plot(ntree_fit)
ntree.err <- data.frame(ntree_fit$mse)
# set ntree = 2000

# -----final model-----
set.seed(119)
train_fit <- randomForest(Age~., data = Age_CAG_train, mtry = (ncol(Age_CAG_train) - 1), ntree = 2000,
                          importance = TRUE, proximity = TRUE)
train.imp <- train_fit$importance

# fit microbial age
pred <- data.frame(Agedata$Age, train_fit$predicted); names(pred) <- c('Age', 'Microbial_Age')
lm.pred <- lm(Microbial_Age ~ Age, data = pred)
summary(lm.pred)
predict <- predict(lm.pred, interval = "confidence")
pred <- cbind(pred, predict)

# train and save male RF model and female RF model separately
if (sex == '0') {
  train_fit_female <- train_fit
  train.imp_female <- train.imp
  pred_female <- pred
} else if (sex == '1') {
  train_fit_male <- train_fit
  train.imp_male <- train.imp
  pred_male <- pred
}

# merge microbial age data
pred_female$Sex <- '0'
pred_male$Sex <- '1'
pred_all <- rbind(pred_female, pred_male)
pred_all <- pred_all[order(rownames(pred_all)), ]


# -----difference of microbial age under disease status-----

# adjustment of disease classification
for (i in 18:33) {
  Pheno[, i] <- as.numeric(Pheno[, i])
}
rm(i)
Pheno$Obesity [Pheno$Obesity == 1] <- 0
Pheno$Obesity [Pheno$Obesity == 2] <- 1
Pheno$AMI_CHF_AF_Stroke [Pheno$AMI == 1 | 
                           Pheno$CHF == 1 |
                           Pheno$AF == 1 |
                           Pheno$Stroke == 1] <- 1
Pheno$Gout_HL_FL [Pheno$Gout == 1 |
                    Pheno$Hyperlipidemia == 1 |
                    Pheno$FattyLiver == 1] <- 1
Pheno$Health [Pheno$Rlt_health == 2] <- 1 

# merge data
Pheno_CAG <- cbind(Pheno, CAG_table); rm(Pheno, CAG_table)

# regression with microbial age within disease group
Coef_female <- data.frame(est = 0, stderr = 0, tvalue = 0, pvalue = 0, Dis = 'Dis', Sex = 'Sex', n = 0,
                          est.dis = 0, stderr.dis = 0, tvalue.dis = 0, pvalue.dis = 0, pvalue.dis.ia = 0)
Coef_male <- data.frame(est = 0, stderr = 0, tvalue = 0, pvalue = 0, Dis = 'Dis', Sex = 'Sex', n = 0,
                          est.dis = 0, stderr.dis = 0, tvalue.dis = 0, pvalue.dis = 0, pvalue.dis.ia = 0)

Dis.n <- c(33,30,18,19,40,31,20,39)

for (n in 1:length(Dis.n)) {
  i = Dis.n[n]
  Dis_logic <- (is.na(Pheno_CAG[,i]) == FALSE & Pheno_CAG[,i] == 1)
  Dis.pop <- Pheno_CAG[Dis_logic, ]
  Dis_name <- names(Pheno_CAG)[i]
  
  # for women
  # fit microbial age
  predict <- predict(train_fit_female, newdata = Dis.pop[Dis.pop$Sex == '0', ])
  Dis.pred.female <- data.frame(Dis.pop[Dis.pop$Sex == '0', c(4,3)], predict)
  names(Dis.pred.female)[3] <- 'Gut_Age'; rm(predict)
  # regression
  lm.pred <- lm(Gut_Age ~ Age, data = Dis.pred.female)
  coef0 <- data.frame(summary(lm.pred)$coefficients)
  names(coef0) <- c('est', 'stderr', 'tvalue', 'pvalue')
  coef0$Dis <- Dis_name
  coef0$Sex <- '0'
  coef0$n <- nrow(Dis.pop[Dis.pop$Sex == '0', ])
  coef0 <- coef0[2,]
  predict <- predict(lm.pred, interval = "confidence")
  Dis.pred.female <- cbind(Dis.pred.female, predict); rm(predict, lm.pred)
  Dis.pred.female$Dis <- '1'
  
  # for men
  # fit microbial age
  predict <- predict(train_fit_male, newdata = Dis.pop[Dis.pop$Sex == '1', ])
  Dis.pred.male <- data.frame(Dis.pop[Dis.pop$Sex == '1', c(4,3)], predict)
  names(Dis.pred.male)[3] <- 'Gut_Age'; rm(predict)
  # regression
  lm.pred <- lm(Gut_Age ~ Age, data = Dis.pred.male)
  coef1 <- data.frame(summary(lm.pred)$coefficients)
  names(coef1) <- c('est', 'stderr', 'tvalue', 'pvalue')
  coef1$Dis <- Dis_name
  coef1$Sex <- '1'
  coef1$n <- nrow(Dis.pop[Dis.pop$Sex == '1', ])
  coef1 <- coef1[2,]
  predict <- predict(lm.pred, interval = "confidence")
  Dis.pred.male <- cbind(Dis.pred.male, predict); rm(predict, lm.pred)
  Dis.pred.male$Dis <- '1'
  
  # merge men and women microbial age
  pred_all$Dis <- '0'
  Dis.pred.all <- rbind(pred_all, Dis.pred.female, Dis.pred.male)
  Dis.pred.all <- within(Dis.pred.all, {
    group <- NA
    group [Sex == '0' & Dis == '0'] <- 'f0'
    group [Sex == '0' & Dis == '1'] <- 'f1'
    group [Sex == '1' & Dis == '0'] <- 'm0'
    group [Sex == '1' & Dis == '1'] <- 'm1'
  })
  
  # compare microbial age between healthy men and men with NCD
  lm.fit.Dis <- lm(Gut_Age ~ Age + Dis, data = Dis.pred.all[Dis.pred.all$Sex == '1', ])
  coef.dis <- data.frame(summary(lm.fit.Dis)$coefficients)
  names(coef.dis) <- c('est.dis', 'stderr.dis', 'tvalue.dis', 'pvalue.dis')
  coef1 <- cbind(coef1, coef.dis[3,])
  # compare microbial age between healthy women and women with NCD
  lm.fit.Dis <- lm(Gut_Age ~ Age + Dis, data = Dis.pred.all[Dis.pred.all$Sex == '0', ])
  coef.dis <- data.frame(summary(lm.fit.Dis)$coefficients)
  names(coef.dis) <- c('est.dis', 'stderr.dis', 'tvalue.dis', 'pvalue.dis')
  coef0 <- cbind(coef0, coef.dis[3,])
  
  # interaction with disease in men
  lm.fit.Disia <- lm(Gut_Age ~ Age + Dis + Age*Dis, data = Dis.pred.all[Dis.pred.all$Sex == '1', ])
  coef0$pvalue.dis.ia <- summary(lm.fit.Disia)$coefficients[4,4]
  # interaction with disease in women
  lm.fit.Disia <- lm(Gut_Age ~ Age + Dis + Age*Dis, data = Dis.pred.all[Dis.pred.all$Sex == '0', ])
  coef1$pvalue.dis.ia <- summary(lm.fit.Disia)$coefficients[4,4]
  
  # merge all Coef data
  Coef_female <- rbind(Coef_female, coef0)
  Coef_male <- rbind(Coef_male, coef1)
}

Coef_female <- Coef_female[-1,]
Coef_male <- Coef_male[-1,]
Coef <- rbind(Coef_female, Coef_male)

