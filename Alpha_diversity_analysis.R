library(openxlsx)

# Phenotype data
Pheno <- read.xlsx('Pheno.xlsx', rowNames = TRUE)
Pheno$Sex <- as.factor(Pheno$Sex)
Pheno$Rlt_health <- as.factor(Pheno$Rlt_health)
# Alpha diversity data
  Alpha <- read.xlsx('Alpha_CAG.xlsx', rowNames = TRUE)
  Alpha <- read.xlsx('Alpha_allOTU.xlsx', rowNames = TRUE)
  Alpha <- read.xlsx('Alpha_genus.xlsx', rowNames = TRUE)

# keep covariates and merge data
Agedata <- Pheno[c('Age', 'Sex', 'Rlt_health')]
inmeta <- rownames(Alpha) %in% rownames(Agedata)
Alphadata <- Alpha[inmeta,]; rm(inmeta)
Age_Alpha <- cbind(Agedata, Alphadata)
rm(Agedata, Alphadata)

# keep healthy subjects
Age_Alpha <- Age_Alpha[Age_Alpha$Rlt_health == 1,]

# -----Alpha—Age regression-----
# in healthy men
summary(lm(chao1 ~ Age, data = Age_Alpha[Age_Alpha$Sex == 1,]))
summary(lm(pielou ~ Age, data = Age_Alpha[Age_Alpha$Sex == 1,]))
# in healthy women
summary(lm(chao1 ~ Age, data = Age_Alpha[Age_Alpha$Sex == 0,]))
summary(lm(pielou ~ Age, data = Age_Alpha[Age_Alpha$Sex == 0,]))

# interaction with sex
summary(lm(chao1 ~ Age + Sex + Age*Sex, data = Age_Alpha))
summary(lm(pielou ~ Age + Sex + Age*Sex, data = Age_Alpha))



# -----variance of alpha diversity explained by CAG-----
Alpha_CAG <- read.xlsx('Alpha_CAG.xlsx', rowNames = TRUE)
CAG_table <- read.xlsx('CAG_table.xlsx', rowNames = T)
# CAG abundance clr-transformed
for (i in 1:130) {
  CAG_table[,i] [CAG_table[,i] == 0] <- 0.5
}
rm(i)
CAG_table <- data.frame(compositions::clr(CAG_table))


# merge Alpha-CAG, chao1 | pielou
Alpha_with_CAG <- cbind(Pheno[c('Sex', 'Rlt_health')], Alpha_CAG[c('chao1', 'pielou')], CAG_table)
# select healthy subjects
Alpha_with_CAG <- Alpha_with_CAG[Alpha_with_CAG$Rlt_health == '1', ]

# calculate explained variance
var.expl.CAG <- data.frame(R2 = 0, cor.est = 0, p.cor = 0)
for (i in 1:130) {
  Loop <- Alpha_with_CAG[, c(1:4, i+4)]
  names(Loop)[5] <- 'CAG_i'
  lm.fit <- lm(pielou ~ CAG_i, data = Loop)
  variance <- summary(lm.fit)
  R2 <- variance$r.squared
  cor <- cor.test(Loop$pielou, Loop$CAG_i)
  cor.est <- cor$estimate; p.cor <- cor$p.value
  
  var.expl <- data.frame(R2, cor.est, p.cor)
  rownames(var.expl) <- paste0('CAG_', i)
  var.expl.CAG <- rbind(var.expl.CAG, var.expl)
  rm(Loop, lm.fit, variance, R2, cor, cor.est, p.cor, var.expl)
}
var.expl.CAG <- var.expl.CAG[-1,]
var.expl.CAG$var.expl <- var.expl.CAG$R2 * 100
var.expl.CAG$p.BH <- p.adjust(var.expl.CAG$p.cor, method = 'BH')


# -----difference under disease status-----
# merge data
Dis_Alpha <- cbind(Pheno, Alpha[rownames(Pheno), ])
for (i in 18:33) {
  Dis_Alpha[, i] <- as.numeric(Dis_Alpha[, i])
}
rm(i)
# Adjustment of disease classification
Dis_Alpha$Obesity [Dis_Alpha$Obesity == 1] <- 0
Dis_Alpha$Obesity [Dis_Alpha$Obesity == 2] <- 1
Dis_Alpha$AMI_CHF_AF_Stroke [Dis_Alpha$AMI == 1 | 
                             Dis_Alpha$CHF == 1 |
                             Dis_Alpha$AF == 1 |
                             Dis_Alpha$Stroke == 1] <- 1
Dis_Alpha$Gout_HL_FL [Dis_Alpha$Gout == 1 |
                      Dis_Alpha$Hyperlipidemia == 1 |
                      Dis_Alpha$FattyLiver == 1] <- 1
Dis_Alpha$Health [Dis_Alpha$Rlt_health == 2] <- 1 

# coefficient by sex
Coef <- data.frame(est = 0, stderr = 0, tvalue = 0, pvalue = 0,
                   Dis = 'Dis', Sex = 'Sex', n = 0, alpha = 'alpha')
for (sex in c('0', '1')) {
  for (i in c(47,33,30,18,19,46,31,20,45)) {
    Dis_logic <- is.na(Dis_Alpha[,i]) == FALSE & Dis_Alpha[,i] == 1
    Dis <- Dis_Alpha[Dis_logic & Dis_Alpha$Sex == sex, ]
    coef<- data.frame(summary(lm(chao1 ~ Age, data = Dis))$coefficients)
    names(coef) <- c('est', 'stderr', 'tvalue', 'pvalue')
    coef$Dis <- names(Dis_Alpha)[i]
    coef$Sex <- sex
    coef$n <- nrow(Dis)
    coef$alpha <- 'chao1'
    coef <- coef[2,]
    Coef <- rbind(Coef, coef)
    rm(Dis_logic, Dis, coef)
    
    Dis_logic <- is.na(Dis_Alpha[,i]) == FALSE & Dis_Alpha[,i] == 1
    Dis <- Dis_Alpha[Dis_logic & Dis_Alpha$Sex == sex, ]
    coef<- data.frame(summary(lm(pielou ~ Age, data = Dis))$coefficients)
    names(coef) <- c('est', 'stderr', 'tvalue', 'pvalue')
    coef$Dis <- names(Dis_Alpha)[i]
    coef$Sex <- sex
    coef$n <- nrow(Dis)
    coef$alpha <- 'pielou'
    coef <- coef[2,]
    Coef <- rbind(Coef, coef)
    rm(Dis_logic, Dis, coef)
  }
}
Coef <- Coef[-1,]

# interaction with sex
IA <- data.frame(est = 0, stderr = 0, tvalue = 0, pvalue = 0,
                   Dis = 'Dis', Sex = 'Sex', n = 0, alpha = 'alpha')
for (sex in c('0', '1')) {
  for (i in c(33,30,18,19,46,31,20,45)) {
    Dis_logic <- (is.na(Dis_Alpha[,i]) == FALSE & Dis_Alpha[,i] == 1) | Dis_Alpha[,33] == 2
    Dis <- Dis_Alpha[Dis_logic & Dis_Alpha$Sex == sex, ]
    Dis$Disgroup [Dis[,33] == 2] <- '0'
    Dis$Disgroup [Dis[,33] != 2] <- '1'
    
    ia<- data.frame(summary(lm(chao1 ~ Age + Disgroup + Age*Disgroup, data = Dis))$coefficients)
    names(ia) <- c('est', 'stderr', 'tvalue', 'pvalue')
    ia$Dis <- names(Dis_Alpha)[i]
    ia$Sex <- sex
    ia$n <- nrow(Dis)
    ia$alpha <- 'chao1'
    ia <- ia[4,]
    IA <- rbind(IA, ia)
    rm(Dis_logic, Dis, ia)
    
    Dis_logic <- (is.na(Dis_Alpha[,i]) == FALSE & Dis_Alpha[,i] == 1) | Dis_Alpha[,33] == 2
    Dis <- Dis_Alpha[Dis_logic & Dis_Alpha$Sex == sex, ]
    Dis$Disgroup [Dis[,33] == 2] <- '0'
    Dis$Disgroup [Dis[,33] != 2] <- '1'
    
    ia<- data.frame(summary(lm(pielou ~ Age + Disgroup + Age*Disgroup, data = Dis))$coefficients)
    names(ia) <- c('est', 'stderr', 'tvalue', 'pvalue')
    ia$Dis <- names(Dis_Alpha)[i]
    ia$Sex <- sex
    ia$n <- nrow(Dis)
    ia$alpha <- 'pielou'
    ia <- ia[4,]
    IA <- rbind(IA, ia)
    rm(Dis_logic, Dis, ia)
  }
}
IA <- IA[-1,]




