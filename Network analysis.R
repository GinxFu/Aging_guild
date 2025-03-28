library(openxlsx)

# -----1. Similarity of Guild structures-----
library(aricode)

group <- c('',
           'age50','age60','age70','age80',
           'men',
           'men_age50','men_age60','men_age70','men_age80',
           'women',
           'women_age50','women_age60','women_age70','women_age80',
           'healthy',
           'healthy_age50','healthy_age60','healthy_age70','healthy_age80',
           'unhealthy',
           'unhealthy_age50','unhealthy_age60','unhealthy_age70','unhealthy_age80')

results_NVI <- as.data.frame(matrix(NA, nrow = length(group), ncol = length(group),
                                    dimnames = list(c('All', group[-1]),
                                                    c('All', group[-1]))))
results_NVI$k <- results_NVI$K <- NA
results_NVI$NVI_random_sd <- results_NVI$NVI_random_mean <- results_NVI$NVI_random_max <- NA

for (i in 1:length(group)) {
  for (j in 1:length(group)) {
    S1 <- read.xlsx(paste0('Cluster_',group[i],'.xlsx'))
    K1 <- gsub('k','', names(S1)[2])
    k1 <- length(table(S1[,2]))
    
    S2 <- read.xlsx(paste0('Cluster_',group[j],'.xlsx'))
    K2 <- gsub('k','', names(S2)[2])
    k2 <- length(table(S2[,2]))
    
    NVI <- 1 - NVI(S1[,2], S2[,2])
    
    results_NVI[i,j] <- NVI
    results_NVI$K[i] <- K1
    results_NVI$k[i] <- k1
    
    # 随机分组
    NVI_random <- rep(NA, 1000)
    for (n in 1:1000) {
      set.seed(n)
      NVI_random[n-118] <-
        1 - NVI(S1[,2], 
                S1[,2] [order(rnorm(1477))])
    }
    results_NVI$NVI_random_max[i] <- max(NVI_random)
    results_NVI$NVI_random_mean[i] <- mean(NVI_random)
    results_NVI$NVI_random_sd[i] <- sd(NVI_random)
    
    rm(S1,S2, K1,K2, k1,k2, NVI)
    rm(NVI_random)
  }
  
}

write.xlsx(results_NVI, 'NVI_by_groups.xlsx',rowNames = T)


# -----2. Network density-----
library(igraph)

# read clustering information of 130 guilds
# read SparCC correlation data
Cluster_130guilds <- read.xlsx('Cluster_130guilds.xlsx')
cor_SparCC <- read.csv('cor_SparCC.csv', row.names = 1)


# weighted adjacent matrix
cor_weighted <- cor_SparCC
for (i in 1:ncol(cor_weighted)) {
  cor_weighted[,i] [abs(cor_weighted[,i]) < 0.4] <- 0
}

Edges <- (sum(cor_weighted != 0) - 1477)/2     ## total edges
Weighted.Edges <- (sum(abs(cor_weighted)) - 1477)/2     ## total weighted edges


# inter- and intra- guild edges
Guild_ID <- unique(Cluster_130guilds$Guild)

Intra.Edges <- rep(NA, 130)
Inter.Edges <- rep(NA, 130)

for (i in 1:length(Guild_ID)) {
  OTU_ID <- Cluster_130guilds$OTU_ID [Guild$Guild == Guild_ID[i]]
  cor_edge_loop <- cor_weighted[OTU_ID, OTU_ID]
  sum <- (sum(abs(cor_edge_loop)) - length(OTU_ID))/2
  Intra.Edges <- Intra.Edges + sum     ## intra-guild edges
  
  rm(OTU_ID, cor_edge_loop, sum)
}

Inter.Edges <- Weighted.Edges - Intra.Edges     ## inter-guild edges
