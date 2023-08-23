library(openxlsx)
library(vegan)

#----------import OTU_table or CAG_table or genus_table----------
OTU_table <- read.xlsx('OTU_table.xlsx', rowNames = T)
CAG_table <- read.xlsx('CAG_table.xlsx', rowNames = T)
genus_table <- read.xlsx('genus_table.xlsx', rowNames = T)


# choose one level as unit
Unit_table <- OTU_table
# or
Unit_table <- CAG_table
# or
Unit_table <- CAG_table

#----------Flattening----------
colSums(Unit_table) 
Unit_Flattening <- as.data.frame(
  t(
    rrarefy(t(Unit_table), min(colSums(Unit_table)))
  )
)
colSums(Unit_Flattening) 

Unit <- t(Unit_Flattening)
#----------alpha----------
#richness
richness <- as.data.frame(t(estimateR(Unit)))
names(richness)[1:5] <- c('obs', 'chao1', 'chao1_se', 'ACE', 'ACE_se')
#shannon
shannon <- as.data.frame(diversity(Unit, index = 'shannon', base = 2))
names(shannon)[1] <- 'shannon'
#Pielou
  #SR
SR <- specnumber(Unit, MARGIN = 1)
Pielou <- shannon/SR
names(Pielou)[1] <- 'pielou'
#goods_coverage
goods_coverage <- as.data.frame(1 - rowSums(Unit == 1) / rowSums(Unit))
names(goods_coverage)[1] <- 'goods_coverage'
# merge
Alpha <- data.frame(richness, shannon, Pielou, goods_coverage)[,c(1:2,4,6:8)]
# export
write.xlsx(Alpha, 'Alpha.xlsx', rowNames = T)

#----------beta----------
Bray <- vegdist(Unit, method = 'bray')
Bray_matrix <- as.matrix(Bray)
# export
write.csv(Bray_matrix, 'Bray.csv', row.names = T)


#-------Distance to centroid------
# import data
Bray <- read.csv('Bray_CAG.csv')
rownames(Bray) <- Bray$X
Bray <- Bray[,-1]
Pheno <- read.xlsx('R:\\文档\\6.Research\\2302_Age_CAG\\Pheno.xlsx', rowNames = TRUE)
Pheno$Sex <- as.factor(Pheno$Sex)
Pheno$Rlt_health <- as.factor(Pheno$Rlt_health)

# convert matrix into distance
Bray_dist <- as.dist(Bray)

# calculate distance to centroid
group <- rep('1', nrow(Pheno))
mod <- betadisper(Bray_dist, group, type = c('centroid'))
Bray_DtoC <- data.frame(mod$distances)