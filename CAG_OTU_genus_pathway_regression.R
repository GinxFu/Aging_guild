library(openxlsx)

# import Phenotype data
Pheno <- read.xlsx('Pheno.xlsx', rowNames = T)
Pheno$Sex <- as.factor(Pheno$Sex)
Pheno$Rlt_health <- as.factor(Pheno$Rlt_health)

Agedata <- Pheno[c('Age', 'Sex', 'Rlt_health')]

# import CAG or OTU or genus or pathway data
CAG_table <- read.xlsx('CAG_table.xlsx', rowNames = T)
OTU_table <- read.xlsx('OTU_table.xlsx', rowNames = T)
genus_table <- read.xlsx('genus_table.xlsx', rowNames = T)
pathway_table <- read.xlsx('pathway_table.xlsx', rowNames = T)

# choose one table as analysis table
analysis_table <- CAG_table
# or
analysis_table <- OTU_table
# or
analysis_table <- genus_table
# or
analysis_table <- pathway_table

# clr-transformed
for (i in 1:130) {
  analysis_table[,i] [analysis_table[,i] == 0] <- 0.5
}
rm(i)
analysis_table <- data.frame(compositions::clr(analysis_table))


# keep healthy subjects and merge data
Agedata <- Agedata[Agedata$Rlt_health == '1',]
inmeta <- rownames(analysis_table) %in% rownames(Agedata)
analysisdata <- analysis_table[inmeta,]
rm(inmeta)
# 合并数据
Age_analysis <- cbind(Agedata, analysisdata)


# regression coefficients by sex
# in men
Coef_men <- data.frame(Est = 'Est', Std = 0, t.value = 0, P.value = 0)
for (i in 1:ncol(analysisdata)) {
  Loop <- Age_analysis[Age_analysis$Sex == 1, c(1,i+3)]
  names(Loop)[2] <- 'CAG_i'
  lm.loop <- lm(CAG_i ~ Age, data = Loop)
  result <- summary(lm.loop)$coefficients
  colnames(result) <- c('Est', 'Std', 't.value', 'P.value')
  
  Coef_men <- rbind(Coef_men, result[2,])
  rownames(Coef_men)[i+1] <- colnames(Age_analysis)[i+3]
  rm(i, Loop, lm.loop, result)
}
Coef_men <- Coef_men[-1,]
Coef_men$P.value_BH <- p.adjust(Coef_men$P.value, method = 'BH')

# in women
Coef_women <- data.frame(Est = 'Est', Std = 0, t.value = 0, P.value = 0)
for (i in 1:ncol(analysisdata)) {
  Loop <- Age_analysis[Age_analysis$Sex == 0, c(1,i+3)]
  names(Loop)[2] <- 'CAG_i'
  lm.loop <- lm(CAG_i ~ Age, data = Loop)
  result <- summary(lm.loop)$coefficients
  colnames(result) <- c('Est', 'Std', 't.value', 'P.value')
  
  Coef_women <- rbind(Coef_women, result[2,])
  rownames(Coef_women)[i+1] <- colnames(Age_analysis)[i+3]
  rm(i, Loop, lm.loop, result)
}
Coef_women <- Coef_women[-1,]
Coef_women$P.value_BH <- p.adjust(Coef_women$P.value, method = 'BH')

#-----------interaction with sex--------------------------
IA.pvalue <- data.frame(Unit = 'Unit', pvalue = 0)
for (i in 1:ncol(analysisdata)) {
  Loop <- Age_analysis[, c(1,2,i+3)]
  names(Loop)[3] <- 'CAG_i'
  lm.loop <- lm(CAG_i ~ Age + Sex + Age*Sex, data = Loop)
  p <- data.frame(Unit = colnames(Age_analysis)[i+3], pvalue = anova(lm.loop)[3,5])
  IA.pvalue <- rbind(IA.pvalue, p)
  rm(i, Loop, lm.loop, p)
}
IA.pvalue <- IA.pvalue[-1, ]; rownames(IA.pvalue) <- NULL
IA.pvalue$pvalue_BH <- p.adjust(IA.pvalue$pvalue, method = 'BH')



