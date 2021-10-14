library(data.table)
library(reshape2)
library(dplyr)
library(ggdendro)
library(ggplot2)
library(tibble)
library(xlsx)

#load data ------------------------------------------

ARGsOAP_mapping_df <- fread("E:/ARGsOAP/SARG_mapping_file/SARG_mapping.out")

ARGsOAP.relAbund_df <- fread("E:/NutSync/papers/COPD/Exploratory Data Analysis/6.COPD_vs_Healthy/1.ARGsOAP_relAbund_wideDf.csv")

tmp <- ARGsOAP.relAbund_df %>% reshape2::melt( variable.name="ARGsubtype", value.name = "relAbund") %>% 
  mutate(ARGtype = sapply(ARGsubtype, function(x) ARGsOAP_mapping_df$ARG_type[which(ARGsOAP_mapping_df$ARG_subtype == x)[1] ])) %>%
  group_by(sample, ARGtype) %>%
  summarise(relAbund= sum(relAbund))


ARGtype_rankDf <- tmp %>% group_by(ARGtype) %>% summarise(avgRelAbund = mean(relAbund)) %>% arrange(desc(avgRelAbund))
ARGtypes_ARGsOAP <- c(ARGtype_rankDf$ARGtype[1:10],"other")
otherARGtypes <- ARGtype_rankDf$ARGtype[11:nrow(ARGtype_rankDf)]

tmp$ARGtype[which(tmp$ARGtype %in% otherARGtypes)] <- "other"

ARGtype_relAbund.Df <- tmp %>% group_by(sample, ARGtype) %>% summarise(relAbund = sum(relAbund))
dat.w <- ARGtype_relAbund.Df %>% dcast(sample~ARGtype, value.var = "relAbund") %>% tibble::column_to_rownames("sample")



z <- dat.w


# normalization --------------- 
'normalization is mandatory for cluster analysis'
df <- t(z) 
x <- as.matrix(scale(df))
nor <- t(x)  #每个样品 的 means都将近0， 所有sd都变成1

#clustering --------
distance <- dist(nor)
mydata.hclust<-hclust(distance,method="ward.D")  
plot(mydata.hclust,hang=-1,  labels=rownames(z)) 
dd.row <- as.dendrogram(mydata.hclust)
od <- order.dendrogram(dd.row) 



# grouping--------
member <- cutree(mydata.hclust, 4)




# plotting --------------------


#p2 by ggplots
ddata_x <- dendro_data(dd.row) 
p2<- ggplot(segment(ddata_x)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_dendro() + theme(axis.title.x=element_blank())


# p0 

xx <- scale(df)[, od] 
xx_names <- attr(xx, "dimnames")
xx_names[[2]]


plotDat <- ARGtype_relAbund.Df %>% dplyr::filter(sample %in% rownames(nor))
plotDat$sample <- factor(plotDat$sample,levels = xx_names[[2]])

color_ARGtypes_ARGsOAP.bm_df <- 
  cbind.data.frame(ARGtype = c("multidrug","macrolide-lincosamide-streptogramin","unclassified", "beta-lactam","tetracycline",
                               "bacitracin", "kasugamycin","aminoglycoside","fosmidomycin","vancomycin","other"),
                   color = c("#DC0000FF","#3C5488FF", "#B09C85FF","#D595A7FF","#00A087FF","#F39B7FFF",
                             "#4DBBD5FF","#8491B4FF", "#91D1C2FF","#7E6148FF","gray"),
                   stringsAsFactors=F)

Colors <- sapply(unique(plotDat$ARGtype),function(x) color_ARGtypes_ARGsOAP.bm_df$color[color_ARGtypes_ARGsOAP.bm_df$ARGtype == x])
p0<-ggplot(plotDat, aes(x = sample, y = relAbund, fill = ARGtype)) + 
  geom_col() +  
  scale_fill_manual(values = Colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  xlab("") + ylab("percentage of drug types") + guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 


library(grid)
library(ggpubr)
grid.newpage()
print(p2, vp=viewport(width=1, height=0.2, x=0.5, y=0.87))
print( p0+theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()), 
       vp=viewport(width=0.92, height=0.7, x=0.5, y=0.41))






# p3
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank()
)



sample_clusterNo_df <- cbind.data.frame(
  sample= names(member),
  clusterNo = as.character(member),
  stringsAsFactors=F
)

sample_clusterNo_df$sample <- factor(sample_clusterNo_df$sample, levels = xx_names[[2]])

p3 <- ggplot(sample_clusterNo_df) +
  # geom_point(aes(x=sample, y=0, fill=clusterNo), shape=21, size=3,color="black") +
  geom_text(aes(x=sample, y=0, label=clusterNo, color=clusterNo)) +
  theme_none + theme(legend.position = "none")
p3 



grid.newpage()
print(p3 + theme(legend.position = "none"), vp=viewport(width = 0.92, height = 0.8, x=0.5, y = 0.72 ))
print(p2, vp=viewport(width=1, height=0.2, x=0.5, y=0.87))
print(p0+theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()), 
       vp=viewport(width=0.92, height=0.7, x=0.5, y=0.36))





