library(openxlsx)
library(vegan)

# import Bray data
Bray <- read.csv('Bray_CAG.csv')
Bray <- read.csv('Bray_allOTU.csv')
Bray <- read.csv('Bray_genus.csv')
rownames(Bray) <- Bray$X
Bray <- Bray[,-1]
# import phenotype data
Pheno <- read.xlsx('Pheno.xlsx', rowNames = TRUE)
Pheno$Sex <- as.factor(Pheno$Sex)
Pheno$Rlt_health <- as.factor(Pheno$Rlt_health)
Agedata <- Pheno[c('Age', 'Sex', 'Rlt_health')]
# keep healthy subjects
Agedata <- Agedata[Agedata$Rlt_health == 1, ]

# -----PERMANOVA-----

# keep healthy men and merge data
Agedata_male <- Agedata[Agedata$Sex == 1,]
Bray_matrix <- as.matrix(Bray)
inmeta <- rownames(Bray_matrix) %in% rownames(Agedata_male)
Bray_matrix <- Bray_matrix[inmeta, inmeta]
rm(inmeta)
# Adonis/PERMANOVA
check <- rownames(Agedata_male)%in%colnames(Bray_matrix); table(check); rm(check)
set.seed(119)
RD <- adonis2(formula = Bray_matrix ~ Age,
              data = Agedata_male,
              permutaions = 999, method = "bray")
RD
rm(RD)

# keep healthy women and merge data
Agedata_female <- Agedata[Agedata$Sex == 0,]
Bray_matrix <- as.matrix(Bray)
inmeta <- rownames(Bray_matrix) %in% rownames(Agedata_female)
Bray_matrix <- Bray_matrix[inmeta, inmeta]
rm(inmeta)
# Adonis/PERMANOVA
check <- rownames(Agedata_female)%in%colnames(Bray_matrix); table(check); rm(check)
set.seed(119)
RD <- adonis2(formula = Bray_matrix ~ Age,
              data = Agedata_female,
              permutaions = 999, method = "bray")
RD
rm(RD)


# interaction with sex
Bray_matrix <- as.matrix(Bray)
inmeta <- rownames(Bray_matrix) %in% rownames(Agedata)
Bray_matrix <- Bray_matrix[inmeta, inmeta]
rm(inmeta)
# Adonis/PERMANOVA
check <- rownames(Agedata)%in%colnames(Bray_matrix); table(check); rm(check)
set.seed(119)
RD <- adonis2(formula = Bray_matrix ~ Age + Sex + Age*Sex,
              data = Agedata,
              permutaions = 999, method = "bray")
RD
rm(RD)



# -----PCOA-----
library(ggplot2)
library(ggpubr)
PCOA2 <- cmdscale(Bray_matrix, eig = T, k = 25) # Bray_matrix of all healthy subjects
# Axis
round(PCOA2$eig[1] / sum(PCOA2$eig), 2) * 100
round(PCOA2$eig[2] / sum(PCOA2$eig), 2) * 100
sum(PCOA2$eig[1:25]) / sum(PCOA2$eig)
# score
pcoa2_point <- as.data.frame(PCOA2$points)
pcoa2_eig <- PCOA2$eig
# merge agedata
pcoa2_data <- cbind(Agedata, pcoa2_point)
rm(pcoa2_point)

#-----PERMANOVA variance explained by CAG--------------------
CAG_table <- read.xlsx('CAG_table.xlsx', rowNames = T)
# clr-transformed
for (i in 1:130) {
  CAG_table[,i] [CAG_table[,i] == 0] <- 0.5
}
rm(i)
CAG_table <- data.frame(compositions::clr(CAG_table))
# keep healthy subjects
CAG_table <- CAG_table[rownames(CAG_table) %in% rownames(Agedata),]

#-----Adonis/PERM-ANOVA-----
check <- rownames(CAG_table)%in%colnames(Bray_matrix); table(check); rm(check)

Bray_var <- data.frame(var.expl = 0, pvalue = 0)
for (i in 1:130) {
  Loop <- CAG_table
  names(Loop)[i] <- 'CAG_i'
  set.seed(119)
  RD <- adonis2(formula = Bray_matrix ~ CAG_i,
                data = Loop, permutaions = 999, method = "bray")
  var.expl <- 100*(RD[1,3] / RD[3,3])
  pvalue <- RD[1,5]
  var <- data.frame(var.expl = var.expl, pvalue = pvalue)
  rownames(var) <- paste0('CAG_',i)
  Bray_var <- rbind(Bray_var, var)
  rm(Loop, RD, var.expl, pvalue, var)
}
Bray_var <- Bray_var[-1,]

#----difference under disease status-----
# merge data
Dis_Pheno <- Pheno
for (i in 18:33) {
  Dis_Pheno[, i] <- as.numeric(Dis_Pheno[, i])
}
rm(i)
# adjustment of disease classification
Dis_Pheno$Obesity [Dis_Pheno$Obesity == 1] <- 0
Dis_Pheno$Obesity [Dis_Pheno$Obesity == 2] <- 1
Dis_Pheno$AMI_CHF_AF_Stroke [Dis_Pheno$AMI == 1 | 
                               Dis_Pheno$CHF == 1 |
                               Dis_Pheno$AF == 1 |
                               Dis_Pheno$Stroke == 1] <- 1
Dis_Pheno$Gout_HL_FL [Dis_Pheno$Gout == 1 |
                        Dis_Pheno$Hyperlipidemia == 1 |
                        Dis_Pheno$FattyLiver == 1] <- 1
Dis_Pheno$Health [Dis_Pheno$Rlt_health == 2] <- 1 

#-----PERMANOVA within disease group-----
RD_all <- data.frame(pvalue = 0, Dis = 'Dis', Sex = 'Sex', n = 0)
for (sex in c('0', '1')) {
  for (i in c(18:20, 30, 31, 33, 39:41)) {
    Dis_logic <- is.na(Dis_Pheno[,i]) == FALSE & Dis_Pheno[,i] == 1
    Dis <- Dis_Pheno[Dis_logic & Dis_Pheno$Sex == sex, ]
    
    Bray_matrix <- as.matrix(Bray)
    inmeta <- rownames(Bray_matrix) %in% rownames(Dis)
    Bray_matrix <- Bray_matrix[inmeta, inmeta]
    rm(inmeta)
    check <- rownames(Dis)%in%colnames(Bray_matrix); table(check); rm(check)
    set.seed(119)
    RD <- adonis2(formula = Bray_matrix ~ Age, 
                  data = Dis,
                  permutaions = 999, method = "bray")
    RD$Dis <- names(Dis_Pheno)[i]
    RD$Sex <- sex
    RD$n <- nrow(Dis)
    RD <- RD[1, c(-1:-4)]; names(RD)[1] <- 'pvalue'
    RD_all <- rbind(RD_all, RD)
    rm(Dis_logic, Dis, Bray_matrix, RD)
  }
}
RD_all <- RD_all[-1,]

#-----PERMANOVA interaction with disease-----
RD_ia <- data.frame(pvalue = 0, Dis = 'Dis', Sex = 'Sex')
for (sex in c('0', '1')) {
  for (i in c(18:20, 30, 31, 33, 39,40)) {
    Dis_logic <- (is.na(Dis_Pheno[,i]) == FALSE & Dis_Pheno[,i] == 1) | Dis_Pheno[,33] == 2
    Dis <- Dis_Pheno[Dis_logic & Dis_Pheno$Sex == sex, ]
    Dis$Disgroup [Dis[,33] == 2] <- '0'
    Dis$Disgroup [Dis[,33] != 2] <- '1'
    
    Bray_matrix <- as.matrix(Bray)
    inmeta <- rownames(Bray_matrix) %in% rownames(Dis)
    Bray_matrix <- Bray_matrix[inmeta, inmeta]
    rm(inmeta)
    check <- rownames(Dis)%in%colnames(Bray_matrix); table(check); rm(check)
    set.seed(119)
    RD <- adonis2(formula = Bray_matrix ~ Age + Disgroup + Age*Disgroup, 
                  data = Dis,
                  permutaions = 999, method = "bray")
    RD$Dis <- names(Dis_Pheno)[i]
    RD$Sex <- sex
    RD <- RD[3, c(-1:-4)]; names(RD)[1] <- 'pvalue'
    RD_ia <- rbind(RD_ia, RD)
    rm(Dis_logic, Dis, Bray_matrix, RD)
  }
}
RD_ia <- RD_ia[-1,]

