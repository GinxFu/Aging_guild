
library(SpiecEasi)
library(openxlsx)
library(vegan)

# -----import OTU table | n = 2944 | 1477 OTU with prevalence > 20%-----
OTU_table <- read.xlsx('OTU_table_2944n_1477OTU.xlsx', rowNames = T)
OTU_table_matrix <- as.matrix(OTU_table)
OTU_id <- colnames(OTU_table)

# -----SpiecEasi-sparcc-----
SparCC_model <- sparcc(OTU_matrix_count, iter = 20, inner_iter = 10, th = 0.1)

cor <- SparCC_model$Cor
cor <- as.data.frame(cor)

write.csv(cor, file = 'cor.csv')

# -----Ward clustering-----
cor <- read.csv('cor.csv')
# convert into distance metrics
cor[c(2:1478)] <- 1 - cor[c(2:1478)]
cor <- cor[,c(2:1478)]
cor_dist <- as.dist(cor)

wardhc <- hclust(cor_dist, method = "ward.D")

# -----PERMANOVA test clustering tree-----

# clustering grouping information 
K_f = function(n){
  K <- data.frame(V = c(1:1477))
  for (i in 1:n) {
    k <- as.data.frame(cutree(wardhc, i))
    names(k)[1] <- paste('k',i,sep = '')
    K <- cbind(K,k)
  }
  K
}
K_total <- K_f(200)
rownames(K_total) <- OTU_id

# PERMANOVA test as the clustering height decreased
K_dupl_PERMANOVA = function(x){
  K_n <- data.frame(table(K_total[,c(x+1,x+2)]))
  names(K_n)[1] <- 'kx'
  names(K_n)[2] <- 'kx+1'
  K_n <- K_n[K_n$Freq > 0,]
  
  dupl <- data.frame(table(K_n[,1]))
  dupl <- dupl[dupl$Freq == 2,]
  K_n <- K_n[K_n$kx == dupl$Var1,]
  
  wardhc_k <- as.data.frame(cutree(wardhc, x+1))
  names(wardhc_k)[1] <- 'k'
  wardhc_k$k <- as.factor(wardhc_k$k)
  m <- K_n$`kx+1`[1]
  n <- K_n$`kx+1`[2]
  wardhc_dupl <- wardhc_k; wardhc_dupl$V <- rownames(wardhc_k)
  wardhc_dupl <- wardhc_dupl[wardhc_dupl$k == m | wardhc_dupl$k == n,]
  
  rownames(cor) <- colnames(cor)
  t <- wardhc_dupl$V
  cor_dupl <- cor[t,t]
  cor_dupl_matrix <- as.matrix(cor_dupl)
  
  check <- rownames(wardhc_dupl)%in%colnames(cor_dupl_matrix);table(check)
  set.seed(119)
  PERMANOVA_p <- adonis2(formula = cor_dupl_matrix ~ k,
                         data = wardhc_dupl,
                         permutaions = 9999)
    
  PERMANOVA_p
}

P <- data.frame(k_p = 1)
for (j in 1:199) {
  P <- cbind(P, data.frame(K_dupl_PERMANOVA(j)[1,5]))
}
P_t <- data.frame(t(P))
rownames(P_t) <- paste0('k',c(1:200))

rownames(P_t)[P_t$t.P. > 0.05]
# "k1"   "k131" "k148" "k149" "k156" "k159" "k162" "k165" "k183" "k184" "k197"
# As the clustering height decreased, the No.131 clade was not significantly different.
# Thus, 130 CAGs were confirmed.

K_130 <- data.frame(OTU_id = rownames(K_total), CAG = K_total[,'k130'])

