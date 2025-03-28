library(openxlsx)
library(vroom)
library(SpiecEasi)
library(vegan)

# read abundance data of 1477 OTUs
OTU_table_count <- vroom('OTU_table.csv', show_col_types = F)

# -----SparCC correlation-----
set.seed(123)
model_SparCC <- sparcc(data_OTU_table, iter = 20, inner_iter = 10, th = 0.1) 
cor_SparCC <- model_SparCC$Cor
cov_SparCC <- model_SparCC$Cov
dim(cor_SparCC)

# 添加OTU_id
colnames(cor_SparCC) <- colnames(cov_SparCC) <- rownames(cor_SparCC) <- rownames(cov_SparCC) <- 
  rownames(OTU_table_p20_for_Net)

# 导出结果
write.csv(cor_SparCC, file = 'cor_SparCC.csv')


# -----Guild (CAG) construction-----

# distance matrix
cor_SparCC_dist <- 1 - cor_SparCC
cor_SparCC_dist <- as.dist(cor_SparCC_dist)

# Ward clustering
wardhc <- hclust(cor_SparCC_dist, method = "ward.D")

## -----testing-----

# grouping of 200 guilds
K_group = function(n){
  K <- data.frame(OTU = rownames(cor_SparCC))
  for (i in 2:n) {
    k <- as.data.frame(cutree(wardhc, i))
    names(k)[1] <- paste('k',i,sep = '')
    K <- cbind(K,k)
  }
  K
}
K_total <- K_group(200)


# test by decreasing clustering height
cor_SparCC_dist_matrix <- as.matrix(cor_SparCC_dist)

K_dupl_PERMANOVA = function(x){
  K_n <- data.frame(table(K_total[,c(x-1,x)]))
  names(K_n)[1] <- 'kx-1'
  names(K_n)[2] <- 'kx'
  K_n <- K_n[K_n$Freq > 0,]
  
  dupl <- data.frame(table(K_n[,1]))
  dupl <- dupl[dupl$Freq == 2,]
  K_n <- K_n[K_n$`kx-1` == dupl$Var1,]
  
  wardhc_k <- as.data.frame(cutree(wardhc, x))
  names(wardhc_k)[1] <- 'k'
  wardhc_k$k <- as.factor(wardhc_k$k)
  m <- K_n$kx[1]
  n <- K_n$kx[2]
  wardhc_dupl <- wardhc_k; wardhc_dupl$OTU <- rownames(wardhc_dupl)
  wardhc_dupl <- wardhc_dupl[wardhc_dupl$k == m | wardhc_dupl$k == n,]
  
  t <- wardhc_dupl$OTU
  cor_dupl_matrix <- cor_SparCC_dist_matrix[t,t]
  
  check <- rownames(wardhc_dupl)%in%colnames(cor_dupl_matrix);table(check)
  
  set.seed(119)
  PERMANOVA_p <- 
    adonis2(formula = cor_dupl_matrix ~ k,
            data = wardhc_dupl,
            permutaions = 999)
  PERMANOVA_p
}

K_PERMANOVA_pvalue = function(n) {
  P <- data.frame(K = c(1:n), Pvalue = NA)
  for (j in 3:n) {
    P[j,2] <- K_dupl_PERMANOVA(j)[1,5]
  }
  P
} 
P <- K_PERMANOVA_pvalue(200)

# select cut-off point 
for (n in 3:200) {
  if (P[n,2] > 0.05)  {
    K_n <- n-1
    break
  }
}
K_n

K_dupl_PERMANOVA(n)[1,5]
K_dupl_PERMANOVA(K_n)[1,5]

K_use <- K_total[c('OTU',paste0('k',K_n))]
K_use[,2] <- paste0('CAG_', K_use[,2])

write.xlsx(K_use,paste0('Cluster.xlsx'))