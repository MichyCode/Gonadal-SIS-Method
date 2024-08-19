options(error=traceback)

rm(list=ls())

require(monocle3,quietly=TRUE,warn.conflicts=FALSE)
require(dplyr,quietly=TRUE,warn.conflicts=FALSE)
require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)
require(ggpubr,quietly=TRUE,warn.conflicts=FALSE)
require(ggthemes,quietly=TRUE,warn.conflicts=FALSE)
require(gridExtra,quietly=TRUE,warn.conflicts=FALSE)
require(reshape2,quietly=TRUE,warn.conflicts=FALSE)
require(MASS)
#require(rayshader)

source('znrf3.r')

genes <- c("Gng13","Pax8", "Amh", "Fgf9", "Sox9", "Sry", "Fst", "Rspo1", "Cyp17a1", "Sry", "Pax8", "Znrf3", "Rnf43",
           "Gata4", "Wt1", "Sf1", "Pax2", "Top2a", "Shisa3", "Xist", "Arhgap29", "Gadd45g", "Foxl2", "Dppa5a", "Fzd1", 
           "Egfl6", "Kctd14", "Wnt6", "Wnt4", "Star", "Insl3", "Irx3", "Col12a1", "Tcf21", "Wnt5a", "Nr2f2", 'Upk3b', 
           'Arx', 'Lhx9', 'Emx2', 'Amhr2', 'Lrrn4', 'Muc16', 'Cbln1', 'Aldh1a1', 'Nr5a1', 'Upk1b', 'Cxcl12', 
           'Cited2', 'Dmrt1', 'Tbx3', 'Top2a')

genes<-unique(genes)

#genes <- c('Axin2','Lef1')

#genes <- c('Amh', 'Sox9', 'Irx3', 'Fst', 'Foxl2', 'Kctd14', 'Upk3b', 'Tcf21', 'Wnt5a', 'Insl3', 'Cyp17a1', 'Dazl', 'Dppa5a', 'Cdh5') 

stageGeneInteractionModel<-F
stepAIC<-F
cData<-rawData

combi<-cData==rawData
sertoli<-grepl("ertoli",rawData,fixed=TRUE)

algorithm<-algName(stepAIC,stageGeneInteractionModel)

cData<-add(cData,paste0("_",algorithm,ifelse(combi,"_confused",ifelse(sertoli,"_sism","_sisf"))))
type<-atype(cData)

message(paste0("Algorithm: ",algorithm,'\t','Dataset type: ',type))

gdir0<-paste0(gdir,algorithm,"/",type,"/")

dummy<-R.utils::mkdirs(gdir0,mustWork=TRUE)

if(!file.exists(cData)) sstop(paste(cData,"does not exist"),1)

message(paste0("\tLoading data from ",cData))

#cData<-'/scratch/deep-znrf3-4/data/znrf3_obj_all_gens_clusters_labelled_with_confused.rds'

data<-readRDS(cData) 
norm_counts <- normalized_counts(data)

get_norm_counts <- function(cds, norm_counts, gene) {
  gene1 <- grab_en_id(gene, cds)
  g1 <- norm_counts[gene1, ] # extract gene row.
  counts <- as.matrix(g1)  # or as.data.frame(as.matrix(obj))
}


for(gene in genes)
{
  #z<-data.frame(data$sism,data$sisf,data$sismf,data$stage,data$cell_type,data$genotype,data$indivlabel,data$confused,colnames(data),data$coex)
  #colnames(z)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","sample","confused","cell","coex")

  # get expression
  x <- get_norm_counts(data, norm_counts,gene)

  z<-data.frame(data$sism,data$sisf,data$sismf,data$stage,data$cell_type,data$genotype,data$confused,colnames(data),data$indiv_label,data$SgeneX5,data$GgeneX5,(data$SgeneX5>=5)&(data$GgeneX5>=5),data$SgeneX5*data$GgeneX5,x)
  colnames(z)<-c("SCIS","GCIS","CCIS","stage","cell_type","genotype","confused","cell","sample","sx","gx","xx","xx2",'x')

  z$stage<-stager(as.factor(z$stage))
  z$genotype<-as.factor(z$genotype)
  z$genotype<-factor(z$genotype,level=c('xx','xy_hom','xy'))

  # Density plots

  z1<-z[z$cell_type %in% interesting,]

  for(ct in c("all cell types","supporting cell types",interesting))
  {
     message(ct)
   
     if(ct=="all cell types")
     {
       z2<-z
     }
     else if(ct=="supporting cell types")
     {
       z2<-z1
     } else
     {
       z2<-z1[z1$cell_type==ct,]
     }
     
     myplots <- list()
 
     xmax<-max(z1$SCIS)
     xmin<-min(z1$SCIS)
     ymax<-max(z1$GCIS)
     ymin<-min(z1$GCIS)
 
     lim<-min(ymin,xmin)
 
     for(stage in sort(unique(z$stage)))
     {
       message(paste0("\t",stage))
     
       for(genotype in c('xx','xy_hom','xy'))
       {
         z3<-z2[(z2$stage==stage)&(z2$genotype==genotype),]

         #z3[is.na(z3$confused_orig)|(z3$SCIS<0)|(z3$GCIS<0),'confused_orig']=FALSE
          
         message(paste0("\t",genotype,"\t",nrow(z3)))

         variation<-sd(z3$x)
  
         if(nrow(z3)<5)
         {
           if(is.na(variation)||(variation==0))
           {
             p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=x)) +
               my_geom_point_no_alpha +
               scale_colour_gradientn(colours = c(baseColors[1]))+
               xlim(xmin,xmax) +
               ylim(ymin,ymax) +
               coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
               mytheme
           } else {
             p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=x)) +
               my_geom_point_no_alpha +
               gene_gradient+
               xlim(xmin,xmax) +
               ylim(ymin,ymax) +
               coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
               mytheme
            }
          } else
          {
            b<-10^(seq(0,-5, length.out=12))

	    if(is.na(variation)||(variation==0))
            {
              p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=x)) +
                my_geom_point_no_alpha +
                scale_colour_gradientn(colours = c(baseColors[1]))+
                #scale_colour_gradientn(colours = terrain.colors(10))+
                xlim(xmin,xmax) +
                ylim(ymin,ymax) +
                coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
                geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+ #,breaks=b)+
                facet_grid(stage~genotype)+
                mytheme
            } else
            {
              p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=x)) +
                my_geom_point_no_alpha +
                #scale_colour_manual(values=c('dark grey','red'))+
                gene_gradient+
                xlim(xmin,xmax) +
                ylim(ymin,ymax) +
                coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
                geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+ #,breaks=b)+
                facet_grid(stage~genotype)+
                mytheme
	     }
          }
    
          myplots[[length(myplots)+1]]<-p  
        }
      }
    
      fn<-paste0(gdir,"by_gene/","contour_",gene,"_",ct,".png")

      message(paste0("\tWriting to ",fn))
    
      p3<-ggarrange(plotlist = myplots,ncol = 3,nrow = 3)
      p3<-annotate_figure(p3,top=ct)
    
      ggsave(file=fn,p3,width=9,height=9,dpi=dpi)
    }
  }
