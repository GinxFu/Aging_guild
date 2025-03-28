library(randomForest)

# read metadata 
data_metadata_healthy <- read.xlsx('metadata_healthy.xlsx')
data_metadata_unhealthy <- read.xlsx('metadata_unhealthy.xlsx')

ID_men <- data_metadata_healthy$ID [data_metadata_healthy$Sex == 1]
ID_women <- data_metadata_healthy$ID [data_metadata_healthy$Sex == 0]

# read abundance data
data_abundance <- read.csv('abundance_table.csv', row.names = 1)
# clr transformation
for (i in 1:nrow(data_abundance_healthy)) {
  data_abundance_healthy[i,] [data_abundance_healthy[i,] == 0] <- 0.5
}
rm(i)
data_abundance <- data.frame(compositions::clr(data_abundance))

data_abundance_healthy <- data_abundance[data_metadata_healthy$ID, ]
data_abundance_unhealthy <- data_abundance[data_metadata_unhealthy$ID, ]


# Significant microbiota in men or women
Micro_sig <- read.xlsx('Abundance_men.xlsx')
Micro_sig <- read.xlsx('Abundance_women.xlsx')

ID_Micro_sig <- Micro_sig$feature [Micro_sig$pval < 0.05]


#------random forest model----------------

# ID_men or ID women
data_RF_train <- cbind(data_metadata_healthy[ID_men, c('Age'), drop = F], 
                       data_abundance_healthy[ID_men, ID_Micro_sig])

#----- 1. Parameter tuning-----

## 1.1 Feature number
set.seed(123)
rfcv.mod.re5 <- replicate(5, rfcv(data_RF_train[,-1], data_RF_train[,1], cv.fold = 10, step = 0.9), simplify = FALSE)
error.cv <- data.frame(sapply(rfcv.mod.re5, '[[', 'error.cv'))
error.cv$err.mean <- rowMeans(error.cv)
# error plot
matplot(rfcv.mod.re5[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type = "l",
        lwd = c(2, rep(1, ncol(error.cv))), col = 1, lty = 1, log = "x",
        xlab = "Number of Feature", ylab = "CV Error")

# Feature number should > 20

## 1.2 mtry
set.seed(123)
fit.mtry <- tuneRF(data_RF_train[,-1], data_RF_train[,1], mtryStart = 1, ntreeTry = 50, stepFactor = 1)
fit.mtry

# set mtry as default: Feature number / 3 

## 1.3 ntree
ntree.err <- data.frame(ntree = c(1:2000), err1 = NA, err2 = NA, err3 = NA, err4 = NA, err5 = NA)
for (i in c(1:5)) {
  set.seed(123+i)
  ntree_fit <- randomForest(Age~., data = data_RF_train, ntree = 2000)
  result_ntree.err[,i+1] <- ntree_fit$mse
}

plot(ntree.err)

# set ntree as 1000 


#-----2. Running model-----

# 2.1 training
set.seed(123)
model_train <- randomForest(Age~., data = data_RF_train, ntree = 1000, importance = TRUE, proximity = TRUE)
model_train
result_train.imp <- as.data.frame(model_train$importance)
varImpPlot(model_train)

write.xlsx(result_train.imp, 'Imp_train.xlsx', rowNames = T)

# get predicted value as Microbial age
result_train_pred <- data.frame(Age = data_RF_train$Age, Micro_age = model_train$predicted)

# linear regression
model_lm.train_pred <- lm(Micro_age ~ Age, data = result_train_pred)
summary(model_lm.train_pred)
result_train_pred <- cbind(result_train_pred, 
                           predict(model_lm.train_pred, interval = "confidence"))


# 2.2 applying
data_RF_apply <- cbind(data_metadata_unhealthy, data_abundance_unhealthy)

result_apply_pred <- predict(model_train, newdata = data_RF_apply)
result_apply_pred <- data.frame(Age = data_RF_apply$Age, Micro_age = result_apply_pred)

# linear
result_apply_pred <- cbind(result_apply_pred,
                          predict(model_lm.train_pred, newdata = result_apply_pred, interval = "confidence"))

# 2.3 merge and output results
result_Micro_age <- rbind(result_train_pred, result_apply_pred)
result_Micro_age <- result_Micro_age[order(rownames(result_Micro_age)), ]
result_Micro_age$ID <- rownames(result_Micro_age)
result_Micro_age <- result_Micro_age[names(result_Micro_age)[c(6,1:5)]]

write.xlsx(result_Micro_age, 'Microbial age men.xlsx', rowNames = F)
write.xlsx(result_Micro_age, 'Microbial age women.xlsx', rowNames = F)
