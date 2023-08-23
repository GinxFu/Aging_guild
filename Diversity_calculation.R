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
