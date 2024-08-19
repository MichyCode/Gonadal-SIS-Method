rm(list=ls())

options(error=traceback)

require(monocle3,quietly=TRUE,warn.conflicts=FALSE)
require(dplyr,quietly=TRUE,warn.conflicts=FALSE)
require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)
require(ggthemes,quietly=TRUE,warn.conflicts=FALSE)
require(gridExtra,quietly=TRUE,warn.conflicts=FALSE)
require(reshape2,quietly=TRUE,warn.conflicts=FALSE)

source('znrf3.r')

cols<-c('green','red','black','blue','blue')
rgbPalette <- c("red","blue","green")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
zPalette <- c('green','red','blue','black')

stage<-NA # 'E12.5'

for(stageGeneInteractionModel in c(FALSE))
{
  for(stepAIC in c(FALSE))
  {
    for(cData in c(rawData)) #,rawData2))
    {
      combi<-cData==rawData
      sertoli<-grepl("ertoli",rawData,fixed=TRUE)
      
      algorithm<-algName(stepAIC,stageGeneInteractionModel)
      
      cData<-add(cData,paste0("_",algorithm,ifelse(combi,"_confused",ifelse(sertoli,"_sism","_sisf"))))
      type<-atype(cData)
    
      message(paste0("Algorithm: ",algorithm,'\t','Dataset type: ',type))
    
      gdir0<-gdir
      
      dummy<-R.utils::mkdirs(gdir0,mustWork=TRUE)
      
      if(!file.exists(cData)) sstop(paste(cData,"does not exist"),1)
        
      message(paste0("\tLoading data from ",cData))

      #cData<-'/scratch/deep-znrf3-4/data/znrf3_obj_all_gens_clusters_labelled_with_confused.rds'
      
      data<-readRDS(cData) 
      
      #z<-data.frame(data$sism,data$sisf,data$sismf,data$stage,data$cell_type,data$genotype,data$indivlabel,data$confused,colnames(data),data$coex)
      #colnames(z)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","sample","confused","cell","coex")

      z<-data.frame(data$sism,data$sisf,data$sismf,data$stage,data$cell_type,data$genotype,data$confused,colnames(data))
      colnames(z)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","confused","cell")
      
      z$stage<-stager(as.factor(z$stage))
      z$genotype<-as.factor(z$genotype)
      z$genotype<-factor(z$genotype,level=c('xx','xy_hom','xy'))
      
      # Density of interesting cell types
      
      message("densities")
      
      z4<-z[z$cell_type %in% interesting,]
      
      youLegend(file=paste0(gdir0,"density_by_cell_type.png"),plot=ggplot(data=z4[z4$genotype=='xy_hom',], aes(x=sis_f, group=cell_type, fill=cell_type)) +
        geom_density(adjust=1.5,alpha=0.4)+ggtitle("Density of sisf by cell type"))

      youLegend(file=paste0(gdir0,"density_by_cell_type_xx.png"),plot=ggplot(data=z4[z4$genotype=='xx',], aes(x=sis_f, group=cell_type, fill=cell_type)) +
               geom_density(adjust=1.5,alpha=0.4)+ggtitle("Density of sisf by cell type - xx"))
      
      youLegend(file=paste0(gdir0,"density_by_cell_type_xy.png"),plot=ggplot(data=z4[z4$genotype=='xy',], aes(x=sis_f, group=cell_type, fill=cell_type)) +
               geom_density(adjust=1.5,alpha=0.4)+ggtitle("Density of sisf by cell type - xy"))
      
      youLegend(file=paste0(gdir0,"density_by_cell_type_xy_hom.png"),plot=ggplot(data=z4[z4$genotype=='xy_hom',], aes(x=sis_f, group=cell_type, fill=cell_type)) +
        geom_density(adjust=1.5,alpha=0.4)+ggtitle("Density of sisf by cell type - xy hom"))


      #
      
      if(!is.na(stage)) z<-z[z$stage==stage,]
      
      p1<-ggplot(data=z,aes(x=sis_m,y=sis_f))+geom_point()+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      p2<-ggplot(data=z,aes(x=sis_m,y=sis_f))+geom_point(aes(color=genotype))+scale_colour_manual(values=rgbPalette)+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      p3<-ggplot(data=z,aes(x=sis_m,y=sis_f))+geom_point(aes(color=stage),size=0.5)+scale_colour_manual(values=zPalette)+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
  
      youLegend(file=paste0(gdir0,"01 sis_f v sis_m all.png"),plot=p1)
      
      message("2")
      
      youLegend(file=paste0(gdir0,"02 sis_f v sis_m all by genotype.png"),plot=p2)
      youLegend(file=paste0(gdir0,"02 sis_f v sis_m all by stage.png"),plot=p3)
      
      p4<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f))+geom_point(aes(color=stage),size=0.5)+scale_colour_manual(values=zPalette)+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))      
      youLegend(file=paste0(gdir0,"02 sis_f v sis_m XX by stage.png"),plot=p4)

      p5<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f))+geom_point(aes(color=stage),size=0.5)+scale_colour_manual(values=zPalette)+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))      
      youLegend(file=paste0(gdir0,"02 sis_f v sis_m XY by stage.png"),plot=p5)

      p6<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f))+geom_point(aes(color=stage),size=0.5)+scale_colour_manual(values=zPalette)+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))      
      youLegend(file=paste0(gdir0,"02 sis_f v sis_m XY hom by stage.png"),plot=p6)
      
      # Tri plot
      
      message("3a")
      
      p1a<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[1])+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      p2a<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[2])+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      p3a<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[3])+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
     
      ggsave(file=paste0(gdir0,"03a sis_f v sis_m all by genotype split",ifelse(is.na(stage),"",paste0("_",stage)),".png"),grid.arrange(p1a,p3a,p2a,ncol=3),width=21,height=7,dpi=dpi)
      
      p1<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f))+geom_point(aes(colour=z[z$genotype=='xx','coex']))+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.position = "none")+scale_color_gradient2(low="green", mid='green', high="red")
      p2<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f))+geom_point(aes(colour=z[z$genotype=='xy','coex']))+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.position = "none")+scale_color_gradient2(low="green", mid='green', high="red")
      p3<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f))+geom_point(aes(colour=z[z$genotype=='xy_hom','coex']))+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.position = "none")+scale_color_gradient2(low="green", mid='green', high="red")

      ggsave(file=paste0(gdir0,"03a coex all by genotype split",ifelse(is.na(stage),"",paste0("_",stage)),".png"),grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7,dpi=dpi)
   
      # Density plots
      
      #message("density by cell types")
      
      #for(i in c(F,T))   
      #{
      #  z3<-z;
      #  
      #  if(i) z3<-z[z$cell_type %in% interesting,]
      #  
      #  xmin<-min(z3$sis_m)
      #  xmax<-max(z3$sis_m)
      #  ymin<-min(z3$sis_f)
      #  ymax<-max(z3$sis_f)      
      #  
      #  for(lim in c(xmin,0,10,max(z3[z3$genotype=='xx','sis_m'])-5))
      #  {
      #    p<-ggplot(data=z3, aes(x = sis_m, y = sis_f)) +
      #      geom_point(aes(size=1)) +
      #      xlim(lim,xmax) +
      #      ylim(lim, 100) +
      #      coord_cartesian(xlim=c(xmin,xmax),ylim=c(ymin,ymax))+
      #      geom_density_2d_filled(alpha = 0.5) +
      #      geom_density_2d(size = 0.25, colour = "black")+
      #      facet_grid(stage~genotype)
      #    theme(legend.position = "none")            

      #    ggsave(file=paste0(gdir0,"data_4 ","03a ",ifelse(i,"interesting ",""),"density lim=",lim," all by genotype split",ifelse(is.na(stage),"",paste0("_",stage)),".png"),
      #           p,width = 21,height=7)
      #  }
      #}
      
      #for(i in c(F,T))   
      #{
      #  z3<-z
      #  
      #  if(i) z3<-z3[z3$cell_type %in% interesting,]
#
#        xmin<-min(z3$sis_m)
#        xmax<-max(z3$sis_m)
#        ymin<-min(z3$sis_f)
#        ymax<-max(z3$sis_f)      
#        
#        for(s in unique(z$sample))
#        {
#          lim<-xmin
#        
#          p1<-
#            ggplot(data=z3[z3$sample==s,], aes(x = sis_m, y = sis_f)) +
#            geom_point() +
#            #xlim(lim,xmax) +
#            #ylim(10, 100) +
#            coord_cartesian(xlim=c(xmin,xmax),ylim=c(ymin,ymax))+
#            geom_density_2d_filled(alpha = 0.5) +
#            geom_density_2d(size = 0.25, colour = "black")+
#            theme(legend.position = "none")
#        
#            ggsave(file=paste0(gdir0,"03a ",s," ",ifelse(i,"interesting ","")," density lim=",lim," all by genotype split",ifelse(is.na(stage),"",paste0("_",stage)),".png"),
#               p1,width =14,height=14)
#        }
#      }
#      
      #stop("!")
    
      message('3c')
  
      #if(!is.na(stage)) stop("Done stage")
     
      # Tri plot with stage shapes
      
      # p1<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[1],aes(shape=stage))+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      # p2<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[2],aes(shape=stage))+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      # p3<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[3],aes(shape=stage))+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      # 
      # ggsave(file=paste0(gdir0,"03c sis_f v sis_m all by genotype split with stages.png"),grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7)

            
      p1a<-ggplot(data=z[(z$genotype=='xx')&(z$stage=='10.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[1],aes(shape=stage))+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())      
      p1b<-ggplot(data=z[(z$genotype=='xx')&(z$stage=='11.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[1],aes(shape=stage))+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p1c<-ggplot(data=z[(z$genotype=='xx')&(z$stage=='12.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[1],aes(shape=stage))+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p1d<-ggplot(data=z[(z$genotype=='xx')&(z$stage=='14.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[1],aes(shape=stage))+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      
      p2a<-ggplot(data=z[(z$genotype=='xy')&(z$stage=='10.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[2],aes(shape=stage))+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())      
      p2b<-ggplot(data=z[(z$genotype=='xy')&(z$stage=='11.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[2],aes(shape=stage))+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p2c<-ggplot(data=z[(z$genotype=='xy')&(z$stage=='12.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[2],aes(shape=stage))+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p2d<-ggplot(data=z[(z$genotype=='xy')&(z$stage=='14.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[2],aes(shape=stage))+ggtitle("xy")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())

      p3a<-ggplot(data=z[(z$genotype=='xy_hom')&(z$stage=='10.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[3],aes(shape=stage))+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p3b<-ggplot(data=z[(z$genotype=='xy_hom')&(z$stage=='11.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[3],aes(shape=stage))+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p3c<-ggplot(data=z[(z$genotype=='xy_hom')&(z$stage=='12.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[3],aes(shape=stage))+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      p3d<-ggplot(data=z[(z$genotype=='xy_hom')&(z$stage=='14.5 dpc'),],aes(x=sis_m,y=sis_f))+geom_point(colour=rgbPalette[3],aes(shape=stage))+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))+theme(legend.title = element_blank())
      
      ggsave(file=paste0(gdir0,"03c sis_f v sis_m all by genotype split with stages.png"),grid.arrange(p1a,p3a,p2a,
                                                                                                       p1b,p3b,p2b,
                                                                                                       p1c,p3c,p2c,
                                                                                                       p1d,p3d,p2d,ncol=3,nrow=4),width = 21,height=21,dpi=dpi)

      message("3b")
      
      # Tri plot with confused

      p1<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f,color=confused))+
        scale_colour_manual(values=c('red','dark red'))+geom_point()+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      p2<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f,color=confused))+
        scale_colour_manual(values=c('blue','dark blue'))+geom_point()+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      p3<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=confused))+
        scale_colour_manual(values=c('green','dark green'))+
        geom_point()+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      ggsave(file=paste0(gdir0,"03b sis_f v sis_m all by genotype split with confused.png"),grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7,dpi=dpi)
      
      # Tri plot with confused+

      for(stage in unique(z$stage))
      {
        z1<-z[z$stage==stage,]
        
        p1<-ggplot(data=z1[z1$genotype=='xx',],aes(x=sis_m,y=sis_f,color=confused))+
          scale_colour_manual(values=c('black','red'))+geom_point()+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
        p2<-ggplot(data=z1[z1$genotype=='xy',],aes(x=sis_m,y=sis_f,color=confused))+
          scale_colour_manual(values=c('black','red'))+geom_point()+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
        p3<-ggplot(data=z1[z1$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=confused))+
          scale_colour_manual(values=c('black','red'))+
          geom_point()+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
          
        #print(grid.arrange(p1,p3,p2,ncol=3))
        #stop()
          
        #ggsave(file=paste0(gdir0,"03b sis_f v sis_m all by genotype split with confused",mode,".png"),grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7)
        
ggsave(file=paste0("/scratch/deep-znrf3-finalfinal/graphs/confused/",stage,".png"),grid.arrange(p1,p3,p2,ncol=3),width=21,height=7,dpi=dpi)
      }
     
      message("4a")
 
      # Tri plot with cell types
      
      cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      p1<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f,color=cell_type))+
              geom_point()+t+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      p2<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f,color=cell_type))+
              geom_point()+t+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      p3<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=cell_type))+
        geom_point(size=0.5)+t+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      youLegend(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with cell types.png"),plot=p3)
      youLegend(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with cell types.png"),plot=grid.arrange(p1,p3,p2,ncol=3),width = 14,height=7)
      
      #z0<-z[z$cell_type %in% c('Granulosa', 'Gonad Progenitor 3', 'Pre-Granulosa', 'Pre-Steroidogenic 4', 'Pre-Steroidogenic 1', 'Gonad Progenitor 2', 'Pre-Supporting Bmyc', 'Pre-Supporting 1', 'Pre-Steroidogenic 3', 'Pre-Supporting 3', 'Pre-Steroidogenic 5', 'Pre-Steroidogenic 2', "Sertoli", "Pre-Sertoli"),]
      z0<-z[z$cell_type %in% interesting,]

      z0$cell_type<-factor(z0$cell_type,interesting)

      p1<-ggplot(data=z0[z0$genotype=='xx',],aes(x=sis_m,y=sis_f,color=cell_type))+
        geom_point()+cts+t+coord_cartesian(xlim=c(0,100),ylim=c(0,100))

      p2<-ggplot(data=z0[z0$genotype=='xy',],aes(x=sis_m,y=sis_f,color=cell_type))+
        geom_point()+cts+t+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))

      p3<-ggplot(data=z0[z0$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=cell_type))+
        geom_point(size=0.5)+cts+t+coord_cartesian(xlim=c(0,100),ylim=c(0,100))

      library(scales)

      ggsave(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with some cell types.png"),grid.arrange(p1,p2,ncol=2),width = 14,height=7)
      ggsave(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with some cell types2.png"),grid.arrange(p1+theme(legend.position="none"),p2+theme(legend.position="none"),ncol=2),width = 14,height=7)
      ggsave(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with some cell types_legend.png"),get_legend(p2),width = 14,height=7)
      
      #youLegend(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with some cell types.png"),plot=p3)
      #youLegend(file=paste0(gdir0,"04a sis_f v sis_m all by genotype split with some cell types.png"),plot=grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7)
      
      youLegend(file=paste0(gdir0,"04a sis_f v sis_m all by XY-hom with some cell types.png"),plot=grid.arrange(p3+scale_colour_colorblind(),ncol=1),width = 7,height=7)

      # Individual cell type graphs
      
      cell_types<-unique(z$cell_type)
      
      for(ct in cell_types)
      {
        z0<-z[z$cell_type==ct,]
        
        f<-nrow(z0[(z0$genotype=='xx')&z0$confused,])/nrow(z[(z$genotype=='xx')&z$confused,])*100
        p1<-ggplot(data=z0[z0$genotype=='xx',],aes(x=sis_m,y=sis_f,color=confused))+geom_point()+scale_colour_manual(values=c('red','dark red'))+ggtitle(paste0("xx - ",ct," n=",nrow(z0[z0$genotype=='xx',])," ",sprintf("%.2f%%",nrow(z0[z0$genotype=='xx',])/nrow(z[z$genotype=='xx',])*100)," (cell types) | dp ",f," %"))+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        
        f<-nrow(z0[(z0$genotype=='xy')&z0$confused,])/nrow(z[(z$genotype=='xy')&z$confused,])*100
        p2<-ggplot(data=z0[z0$genotype=='xy',],aes(x=sis_m,y=sis_f,color=confused))+geom_point()+scale_colour_manual(values=c('blue','dark blue'))+ggtitle(paste0("xy - ",ct," n=",nrow(z0[z0$genotype=='xy',])," ",sprintf("%.2f%%",nrow(z0[z0$genotype=='xy',])/nrow(z[z$genotype=='xy',])*100)," (cell types) | dp ",f," %"))+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        
        f<-nrow(z0[(z0$genotype=='xy_hom')&z0$confused,])/nrow(z[(z$genotype=='xy_hom')&z$confused,])*100
        p3<-ggplot(data=z0[z0$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=confused))+geom_point()+scale_colour_manual(values=c('green','dark green'))+ggtitle(paste0("xy hom - ",ct," n=",nrow(z0[z0$genotype=='xy_hom',])," ",sprintf("%.2f%%",nrow(z0[z0$genotype=='xy_hom',])/nrow(z[z$genotype=='xy_hom',])*100)," (cell types) | dp ",f," %"))+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        
        youLegend(file=paste0(gdir0,"04b sis_f v sis_m ",ct,".png"),plot=p3)
        youLegend(file=paste0(gdir0,"04b sis_f v sis_m ",ct,".png"),plot=grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7)
      }
      
      # Samples
      
      p1<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f,color=sample))+
        geom_point()+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      p2<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f,color=sample))+
        geom_point()+ggtitle("xx")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      p3<-ggplot(data=z[z$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=sample))+
        geom_point(size=0.5)+ggtitle("xy hom")+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
      
      youLegend(file=paste0(gdir0,"05a sis_f v sis_m all by genotype split with samples.png"),plot=p3,width = 21,height=7)
      ggsave(file=paste0(gdir0,"05a sis_f v sis_m all by genotype split with samples.png"),grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7,dpi=dpi)
      
      samples<-unique(z$sample)
      
      for(s in samples)
      {
        z0<-z[z$sample==s,]
        p1<-ggplot(data=z[z$genotype=='xx',],aes(x=sis_m,y=sis_f,color=confused))+
          geom_point()+
          scale_colour_manual(values=c('red','dark red'))+
          ggtitle(paste0("xx - ",s," n=",-1))+
          coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        
        p2<-ggplot(data=z[z$genotype=='xy',],aes(x=sis_m,y=sis_f,color=confused))+
          geom_point()+
          scale_colour_manual(values=c('blue','dark blue'))+
          ggtitle(paste0("xy - ",s," n=",-1))+ 
          coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        
        p3<-ggplot(data=z0[z0$genotype=='xy_hom',],aes(x=sis_m,y=sis_f,color=confused))+
          geom_point()+
          scale_colour_manual(values=c('green','dark green'))+
          ggtitle(paste0("xy hom - ",s," n=",nrow(z[z$genotype=='xy_hom',])," ",sprintf("%.2f%%",nrow(z0[z0$genotype=='xy_hom',])/nrow(z[z0$genotype=='xy_hom',])*100)))+ coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        
        ggsave(file=paste0(gdir0,"05b sis_f v sis_m ",s,".png"),grid.arrange(p1,p3,p2,ncol=3),width = 21,height=7,dpi=dpi)
      }
    }
  }
}

quit(save="no")

# Stats

# COUNTS OVERALL

ps<-c()
qs<-c()

for(pass in 1:2)
{
  i<-1
  
  for(genotype in c('xx','xy'))
  {
    message(genotype)
  
    for(ct in cell_types)
    {
      hom<-z[(z$genotype=='xy_hom')&(z$cell_type==ct),]
      other<-z[(z$genotype==genotype)&(z$cell_type!=ct),]
  
      ldp<-sum(hom[hom$cell_type==ct,]$confused)
      lt<-nrow(hom[hom$cell_type==ct,])
  
      adp<-sum(other[other$cell_type!=ct,]$confused)
      at<-nrow(other[other$cell_type!=ct,])
  
      #if((lt>0)&&(at>0)&&(ldp/lt>adp/at))
      {
        if(pass==1)
        {
          t<-matrix(c(ldp,lt,adp,at),ncol=2)
          rownames(t)<-c('confused','rest')
          colnames(t)<-c(ct,paste0("not ",ct))
          c<-chisq.test(t)
    
          ps<-c(ps,c$p.value)
        }
        else
        {
          message(paste(genotype,ct,qs[i],(ldp/lt)/(adp/at),ifelse(ldp/lt>adp/at,'+','-'),sep='\t'))
          i<-i+1
        }
      }
    }
  }
  
  if(pass==1) qs<-qvalue::qvalue(ps)$qvalues
}

# COUNTS BY STAGE

ps<-c()
qs<-c()

for(pass in 1:2)
{
  i<-1

  for(stage in c('11.5 dpc','12.5 dpc'))
  {
    for(genotype in c('xx','xy'))
    {
      message(genotype)
      
      for(ct in cell_types)
      {
        hom<-z[(z$stage==stage)&(z$genotype=='xy_hom')&(z$cell_type==ct),]
        other<-z[(z$stage==stage)&(z$genotype==genotype)&(z$cell_type!=ct),]
      
        ldp<-sum(hom[hom$cell_type==ct,]$confused)
        lt<-nrow(hom[hom$cell_type==ct,])
      
        adp<-sum(other[other$cell_type!=ct,]$confused)
        at<-nrow(other[other$cell_type!=ct,])
        
        #if((lt>0)&&(at>0)&&(ldp/lt>adp/at))
        {
          if(pass==1)
          {
            t<-matrix(c(ldp,lt,adp,at),ncol=2)
              rownames(t)<-c('confused','rest')
            colnames(t)<-c(ct,paste0("not ",ct))
            c<-chisq.test(t)
            
            ps<-c(ps,c$p.value)
          }
          else
          {
            message(paste(stage,genotype,ct,qs[i],(ldp/lt)/(adp/at),ifelse(ldp/lt>adp/at,'+','-'),sep='\t'))
            i<-i+1
          }
        }
      }
    }
  }
  
  if(pass==1) qs<-qvalue::qvalue(ps)$qvalues
}

# By Stage/sample

ps<-c()
qs<-c()

for(pass in 1:2)
{
  i<-1
  
  for(stage in c('11.5 dpc','12.5 dpc'))
  {
    for(genotype in c('xx','xy'))
    {
      message(genotype)
      
      for(sample in unique(z$sample))
      {
        message(paste0("\t",sample))
      
        for(ct in unique(z$cell_type))
        {
          hom<-z[(z$sample==sample)&(z$stage==stage)&(z$genotype=='xy_hom')&(z$cell_type==ct),]
          other<-z[(z$stage==stage)&(z$genotype==genotype)&(z$cell_type!=ct),]
          
          if(nrow(hom)>0)
          {
            ldp<-sum(hom[hom$cell_type==ct,]$confused)
            lt<-nrow(hom[hom$cell_type==ct,])
            
            adp<-sum(other[other$cell_type!=ct,]$confused)
            at<-nrow(other[other$cell_type!=ct,])
  
            t<-matrix(c(ldp,lt,adp,at),ncol=2)
            rownames(t)<-c('confused','rest')
            colnames(t)<-c(ct,paste0("not ",ct))
            
            #if((lt>0)&&(at>0)&&(ldp/lt>adp/at))
            if(sum(t)>0)
            {
              if(pass==1)
              {
                c<-chisq.test(t)
                
                  ps<-c(ps,c$p.value)
              }
              else
              {
                message(paste(stage,genotype,sample,ct,qs[i],(ldp/lt)/(adp/at),ifelse(ldp/lt>adp/at,'+','-'),sep='\t'))
                i<-i+1
              }
            }
          }
        }
      }
    }
  }
  
  if(pass==1) qs<-qvalue::qvalue(ps)$qvalues
}















#############################################################

rm(ps,qs,genotype,pass,i,ct)

confused_cells<-z[z$confused,'cell']

confusedf<-cellsetf[cellsetf$cell %in% confused_cells,]
f<-cellsetf[!cellsetf$cell %in% confused_cells,]
confusedm<-cellsetm[cellsetm$cell %in% confused_cells,]
m<-cellsetm[!cellsetm$cell %in% confused_cells,]

message("Gene expression levels confused v not confused:\n")

for(g in 10:ncol(confusedf))
{
  levels<-f[,g]
  clevels<-confusedf[,g]
  
  t<-wilcox.test(levels,clevels, alternative = "two.sided")
  p<-t$p.value
  
  message(paste("f",colnames(confusedf)[g],(mean(clevels)-mean(levels))/sd(c(levels,clevels)),p,sep='\t'))
}

for(g in 10:ncol(confusedm))
{
  levels<-m[,g]
  clevels<-confusedm[,g]
  
  t<-wilcox.test(levels,clevels, alternative = "two.sided")
  p<-t$p.value
  
  message(paste("m",colnames(confusedf)[g],(mean(levels)-mean(clevels))/sd(c(levels,clevels)),p,sep='\t'))
}
#data2<-data[,(data$sism>quantile(z$sism,c(0.75)))&(data$sisf>quantile(z$sisf,c(0.75)))]

#genes<-fData(data2)
#m<-match(d$ensembl,rownames(genes))
#z$gene<-genes[m,'gene_short_name']

#z<-data.frame(matrix(0,nrow=ncol(data2),ncol=0))

#z$sism<-data2$sism
#z$sisf<-data2$sisf
#z$cell_type<-data2$cell_type
#z$genotype<-data2$genotype
#z$label<-data2$label

#p<-ggplot(data.frame(z),aes(x=sisf,y=sism,color=cell_type))+geom_point(size = 1, stroke = 0, shape = 16)
#print(p)

