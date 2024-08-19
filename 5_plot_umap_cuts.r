require(monocle3)
require(ggplot2)
require(ggpubr)

source('znrf3.r')

dir<-paste0(graph,"umap/")

combi<-cData==rawData
sertoli<-grepl("ertoli",rawData,fixed=TRUE)

algorithm<-algName(F,F)

cData<-add(rawData,paste0("_",algorithm,"_confused"))
type<-atype(cData)

stageGeneInteractionModel<-F
stepAIC<-F

for(cData in c(rawData)) #,rawData2))
{
  combi<-cData==rawData
  sertoli<-grepl("ertoli",rawData,fixed=TRUE)
  
  algorithm<-algName(stepAIC,stageGeneInteractionModel)
  
  cData<-add(cData,paste0("_",algorithm,ifelse(combi,"_confused",ifelse(sertoli,"_sism","_sisf"))))
  type<-atype(cData)
  
  if(!exists("d")) d<-readRDS(cData)    
  
  dir1<-paste0(gdir,"test_cutoffs2/")
  dummy<-R.utils::mkdirs(dir1,mustWork=TRUE)
  
  for(score in c("sismf"))    
  {
    dir2<-paste0(dir1,score,"/")
    dummy<-R.utils::mkdirs(dir2,mustWork=TRUE)

    message(paste(score,sep='\t'))

    for(mc in seq(0,100,5))    
    {
      fc<-mc
      #for(fc in seq(0,100,10))
      {
        d1<-d[,(d$sism>=mc)&(d$sism>=fc)]
        
        d2<-d1[,d1$genotype=='xx']
        p1<-plot_cells(d2,color_cells_by='sismf_rank',show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
          theme(legend.position = "none",panel.background = element_rect(fill = "lightgray",
                                                                         colour = "lightgray",
                                                                         size = 0.5, linetype = "solid"))+          
          ggtitle("xx")+ylim(-15,15)+xlim(-15,15)
        
        d2<-d1[,d1$genotype=='xy_hom']
        p2<-plot_cells(d2,color_cells_by='sismf_rank',show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
          theme(legend.position = "none",panel.background = element_rect(fill = "lightgray",
                                                                         colour = "lightgray",
                                                                         size = 0.5, linetype = "solid"))+                    
          ggtitle("xy_hom")+ylim(-15,15)+xlim(-15,15)
        
        d2<-d1[,d1$genotype=='xy']
        p3<-plot_cells(d2,color_cells_by='sismf_rank',show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
          theme(legend.position = "none",panel.background = element_rect(fill = "lightgray",
                                                                         colour = "lightgray",
                                                                         size = 0.5, linetype = "solid"))+          
          ggtitle("xy")+ylim(-15,15)+xlim(-15,15)
        

        
        suppressMessages(ggsave(filename = sprintf("%s_%03d_%03d.png",dir2,mc,fc),
                                annotate_figure(ggarrange(p1,p2,p3,nrow=1),top=paste(score,paste0("- sisf>=",fc),paste0(" & sism>=",mc))),width = 12,height=4))
      }
    }
  }
}
