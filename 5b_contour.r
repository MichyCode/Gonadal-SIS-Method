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

cData<-'/scratch/deep-znrf3-finalfinal/data/all_wt_obj_pc55_k19_Normal_confused.rds'

message(paste0("\tLoading data from ",cData))

data<-readRDS(cData) 

#z<-data.frame(data$sism,data$sisf,data$sismf,data$stage,data$cell_type,data$genotype,data$indivlabel,data$confused,colnames(data),data$coex)
#colnames(z)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","sample","confused","cell","coex")

z<-data.frame(data$sism,data$sisf,data$sismf,data$stage,data$cell_type,data$genotype,data$confused,colnames(data),data$indiv_label,data$SgeneX5,data$GgeneX5,(data$SgeneX5>=5)&(data$GgeneX5>=5),data$SgeneX5*data$GgeneX5)
colnames(z)<-c("SCIS","GCIS","CCIS","stage","cell_type","genotype","confused","cell","sample","sx","gx","xx","xx2")

z$stage<-stager(as.factor(z$stage))
z$genotype<-as.factor(z$genotype)
z$genotype<-factor(z$genotype,level=c('xx','xy_hom','xy'))

z$special<-'Others'
z[z$cell_type %in% c('Pre-Sertoli','Sertoli'),'special']<-'Sertoli&Pre-Sertoli'
z[z$cell_type %in% c('Pre-Supporting_1','Pre-Supoorting_2','Granulosa','SLC'),'special']<-'Pre-Supp,SLC&Granulosa'

z$special<-as.factor(z$special)

# Density plots

z1<-z[z$cell_type %in% interesting,]

ssaveRDS(paste0(dataDir,"b6j_z.rds"),z)

for(dp in c("wt","dp","normal","special"))
{
  for(ct in c("all cell types","supporting cell types",interesting))
  {
    if((dp=='special')&&(!ct %in% c('supporting cell types','all cell types'))) next

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

    message(paste(xmin,xmax,ymin,ymax,sep='\t'))

    if(dp=='wt')
    {
      xmin<- min(-12.4951859968301,-13.9953814276731)
      xmax<- max(112.785509554041,137.251724055012)
      ymin<- min(-20.8809000294161,-43.2181547471086)
      ymax<- max(143.606908823565,127.910096543385)
    }

    lim<-min(ymin,xmin)
  
    for(stage in sort(unique(z$stage)))
    {
      if((stage=='14.5 dpc')&&(dp=='wt')) next

      message(paste0("\t",stage))
      
      for(genotype in c('xx','xy_hom','xy'))
      {
	if((genotype=='xy_hom')&&(dp=='wt')) next

        z3<-z2[(z2$stage==stage)&(z2$genotype==genotype),]
        
        #z3[is.na(z3$confused_orig)|(z3$SCIS<0)|(z3$GCIS<0),'confused_orig']=FALSE
           
        message(paste0("\t",genotype,"\t",nrow(z3)))
    
        if(nrow(z3)<5)
        {
          if(dp=='xx')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=xx)) +
              my_geom_point +
              scale_colour_gradientn(colours = terrain.colors(10))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              mytheme
          }
          else if(dp=='sx')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=sx)) +
              my_geom_point +
              scale_colour_gradientn(colours = terrain.colors(10))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              mytheme
          }
          else if(dp=='gx')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=gx)) +
              my_geom_point +
              scale_colour_gradientn(colours = terrain.colors(10))+              
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              mytheme
          }
          else if(dp=='dp')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=confused)) +
              my_geom_point_no_alpha +
              scale_colour_manual(values=c(dotColor,confusedColor))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              mytheme
          }
          else if(dp=='special')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=special)) +
              my_geom_point +
              scale_colour_manual(values=c('grey','red','blue'))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              mytheme
          }
          else # dp=='normal' or dp=='wt'
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS)) +
              my_geom_point +
              scale_colour_manual(values=c(dotColor,confusedColor))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              mytheme
          }
          
        } else
        {
          b<-10^(seq(0,-4, length.out=10))

          if(dp=='xx')        
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=xx2)) +
              my_geom_point +
              #scale_colour_manual(values=c('dark grey','red'))+
              scale_colour_gradientn(colours = terrain.colors(10))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              #geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,alpha = 0.5)+ #,breaks=b) +
              #scale_fill_brewer(palette="Greens")+
              geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+ #,breaks=b)+
              facet_grid(stage~genotype)+
              mytheme
            
          }
          else if(dp=='sx')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=sx)) +
              my_geom_point +
              #scale_colour_manual(values=c('dark grey','red'))+
              scale_colour_gradientn(colours = terrain.colors(10))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              #geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,alpha = 0.5)+ #,breaks=b) +
              #scale_fill_brewer(palette="Greens")+
              geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+ #,breaks=b)+
              facet_grid(stage~genotype)+
              mytheme
          }
          else if(dp=='gx')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=gx)) +
              my_geom_point +
              scale_colour_gradientn(colours = terrain.colors(10))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              #geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,alpha = 0.5)+ #,breaks=b) +
              #scale_fill_brewer(palette="Greens")+
              geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+ #,breaks=b)+
              facet_grid(stage~genotype)+
              mytheme
          }
          else if(dp=='dp')
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=confused)) +
              my_geom_point_no_alpha +
              scale_colour_manual(values=c(dotColor,confusedColor))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,alpha = 0.5,breaks=b)+ # ,breaks=b
              scale_fill_brewer(palette="Blues",direction=-1)+
              #scale_fill_tron()+
              #scale_fill_discrete()+
              geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+
              #facet_grid(stage~genotype)+
              mytheme
          }
          else if(dp=='special')
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=special)) +
              my_geom_point +
              scale_colour_manual(values=c('grey','red','blue'))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,alpha = 0.5,breaks=b)+ # ,breaks=b
              scale_fill_brewer(palette="Blues",direction=-1)+
              #scale_fill_tron()+
              #scale_fill_discrete()+
              geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black",breaks=b)+
              #facet_grid(stage~genotype)+
              mytheme
          else
          {
            p<-ggplot(data=z3, aes(x = SCIS, y = GCIS)) +
              my_geom_point +
              scale_colour_manual(values=c(dotColor,confusedColor))+
              #scale_size_manual(values=c(dotSize))+
              xlim(xmin,xmax) +
              ylim(ymin,ymax) +
              coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
              geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,alpha = 0.5,breaks=b)+
              scale_fill_brewer(palette="Purples",direction=-1)+
              geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = contourSize, colour = contourColor,breaks=b)+
              #facet_grid(stage~genotype)+
              mytheme
          }
        }
    
        myplots[[length(myplots)+1]]<-p  
      }
    }
    
    fn<-paste0(gdir,"contour_",dp,"_",ct,".png")

    message(paste0("\tWriting to ",fn))
   
    if(dp=='wt')
    {
       p3<-ggarrange(plotlist = myplots,ncol = 2,nrow = 2)
    } else
    {
       p3<-ggarrange(plotlist = myplots,ncol = 3,nrow = 3)
    }

    #p3<-annotate_figure(p3,top=ct)
    
    ggsave(file=fn,plot=p3,width=9,height=9,dpi=dpi)

    #fn<-paste0(gdir,"contour_",dp,"_",ct,"_with_confused_3d.png")
    #render_snapshot(fn,clear = TRUE)    
    
    # if(ct=="supporting cell types")
    # {
    #   z3<-z2[(z2$stage=='E14.5')&(z2$genotype!='xy_hom'),]
    #   
    #   p<-ggplot(data=z3, aes(x = SCIS, y = GCIS,color=genotype)) +
    #     geom_point() +
    #     scale_colour_manual(values=c('black','red'))+
    #     xlim(xmin,xmax) +
    #     ylim(ymin,ymax) +
    #     coord_cartesian(xlim=c(lim,xmax),ylim=c(lim,ymax))+
    #     geom_density_2d_filled(data=z3,aes(x=SCIS,y=GCIS,color=genotype),inherit.aes=F,alpha = 0.25)+ #,breaks=b) +
    #     scale_fill_brewer(palette="Greens")+
    #     geom_density_2d(data=z3,aes(x=SCIS,y=GCIS),inherit.aes=F,size = 0.25, colour = "black")+ #,breaks=b)+
    #     #facet_grid(stage~genotype)+
    #     theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))                  
    # 
    # }
    
  }
}
