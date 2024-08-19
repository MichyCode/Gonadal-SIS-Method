rm(list=ls())

options(error=traceback)

require(monocle3,quietly=TRUE,warn.conflicts=FALSE)
require(dplyr,quietly=TRUE,warn.conflicts=FALSE)
require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)
require(ggthemes,quietly=TRUE,warn.conflicts=FALSE)
require(gridExtra,quietly=TRUE,warn.conflicts=FALSE)
require(reshape2,quietly=TRUE,warn.conflicts=FALSE)

source('znrf3.r')

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

for(pass in c(1,2))
{
  cellset1<-readRDS('/scratch/deep-znrf3-finalfinal/models/sex_gene_all_Normal_all_Sertoli_X_Y.rds')
  cellset2<-readRDS('/scratch/deep-znrf3-finalfinal/models/sex_gene_all_Normal_all_Granulosa_X_Y.rds')

  if(pass==1)
  {
    cellset1<-cellset1[cellset1$genotype %in% c('xx','xy'),]
    cellset2<-cellset2[cellset2$genotype %in% c('xx','xy'),]
  } else
  {
    cellset1<-cellset1[!cellset1$genotype %in% c('xx','xy'),]
    cellset2<-cellset2[!cellset2$genotype %in% c('xx','xy'),]
  }

  cs1<-cellset1[,19:ncol(cellset1)]
  colnames(cs1)<-paste0("S ",colnames(cs1))

  cs2<-cellset2[,19:ncol(cellset2)]
  colnames(cs2)<-paste0("G ",colnames(cs2))

  cs<-cbind(cs1,cs2)


  cs4<-cor(cs)
  cs4<-reorder_cormat(cs4)

  require(reshape2)

  cs5<-melt(cs4)

  p<-ggplot(data=cs5,aes(x=Var1,y=Var2,fill=value))+geom_tile()+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    midpoint = 0, limit = c(-1,1), space = "Lab")

  ggsave(file=paste0(gdir,ifelse(pass==1,"wt","hom"),"_gene_heatmap.png"),p,width=8,height=8,dpi=300)
}
