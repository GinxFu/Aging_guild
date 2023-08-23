
library(openxlsx)

# import phenotype and Bray-Curtis data
Pheno <- read.xlsx('Pheno.xlsx', rowNames = TRUE)
Pheno$Sex <- as.factor(Pheno$Sex)
Pheno$Rlt_health <- as.factor(Pheno$Rlt_health)
Bray <- read.csv('Bray_CAG.csv')
rownames(Bray) <- Bray$X
Bray <- Bray[,-1]

Bray <- Bray[rownames(Bray) %in% rownames(Pheno), 
             rownames(Bray) %in% rownames(Pheno)]

# calculate uniqueness in all subjects (n = 2944)
Uniqueness <- data.frame(Uniqueness = rep(0,ncol(Bray)))
rownames(Uniqueness) <- rownames(Bray)
for (i in 1:ncol(Bray)) {
  Uniqueness[i,] <- min(Bray[-i,i])
}
Uniqueness <- Uniqueness[rownames(Uniqueness) %in% rownames(Pheno), ]

# merge data
Uniqdata <- cbind(Pheno, Uniqueness)

# -----regression-----
# in men
summary(lm(Uniqueness ~ Age, data = Uniqdata[Uniqdata$Rlt_health == '1' & Uniqdata$Sex == '1', ]))
# in women
summary(lm(Uniqueness ~ Age, data = Uniqdata[Uniqdata$Rlt_health == '1' & Uniqdata$Sex == '0', ]))

# interaction
summary(lm(Uniqueness ~ Age + Sex + Age* Sex, data = Uniqdata[Uniqdata$Rlt_health == '1', ]))


# -----difference under diasease status-----
for (i in 18:33) {
  Uniqdata[, i] <- as.numeric(Uniqdata[, i])
}
rm(i)
# adjustment of disease classification
Uniqdata$Obesity [Uniqdata$Obesity == 1] <- 0
Uniqdata$Obesity [Uniqdata$Obesity == 2] <- 1
Uniqdata$AMI_CHF_AF_Stroke [Uniqdata$AMI == 1 | 
                               Uniqdata$CHF == 1 |
                               Uniqdata$AF == 1 |
                               Uniqdata$Stroke == 1] <- 1
Uniqdata$Gout_HL_FL [Uniqdata$Gout == 1 |
                        Uniqdata$Hyperlipidemia == 1 |
                        Uniqdata$FattyLiver == 1] <- 1
Uniqdata$Health [Uniqdata$Rlt_health == 2] <- 1 

# coef within disease group
Coef <- data.frame(est = 0, stderr = 0, tvalue = 0, pvalue = 0,
                   Dis = 'Dis', Sex = 'Sex', n = 0)
for (sex in c('0', '1')) {
  for (i in c(43,33,30,18,19,42,31,20,41)) {
    Dis_logic <- is.na(Uniqdata[,i]) == FALSE & Uniqdata[,i] == 1
    Dis <- Uniqdata[Dis_logic & Uniqdata$Sex == sex, ]
    coef<- data.frame(summary(lm(Uniqueness ~ Age, data = Dis))$coefficients)
    names(coef) <- c('est', 'stderr', 'tvalue', 'pvalue')
    coef$Dis <- names(Uniqdata)[i]
    coef$Sex <- sex
    coef$n <- nrow(Dis)
    coef <- coef[2,]
    Coef <- rbind(Coef, coef)
    rm(Dis_logic, Dis, coef)
  }
}
Coef <- Coef[-1,]

# interaction with disease
IA <- data.frame(est = 0, stderr = 0, tvalue = 0, pvalue = 0,
                   Dis = 'Dis', Sex = 'Sex', n = 0)
for (sex in c('0', '1')) {
  for (i in c(33,30,18,19,42,31,20,41)) {
    Dis_logic <- (is.na(Uniqdata[,i]) == FALSE & Uniqdata[,i] == 1) | Uniqdata[,33] == 2
    Dis <- Uniqdata[Dis_logic & Uniqdata$Sex == sex, ]
    Dis$Disgroup [Dis[,33] == 2] <- '0'
    Dis$Disgroup [Dis[,33] != 2] <- '1'
    
    ia <- data.frame(summary(lm(Uniqueness ~ Age + Disgroup + Age*Disgroup, data = Dis))$coefficients)
    names(ia) <- c('est', 'stderr', 'tvalue', 'pvalue')
    ia$Dis <- names(Uniqdata)[i]
    ia$Sex <- sex
    ia$n <- nrow(Dis)
    ia <- ia[4,]
    IA <- rbind(IA, ia)
    rm(Dis_logic, Dis, ia)
  }
}
IA <- IA[-1,]
