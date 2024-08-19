rm(list=ls())

options(error=traceback)

require(monocle3)
require(ggplot2)

source('znrf3.r')

dir<-paste0(graph,"umap/")

combi<-cData==rawData
sertoli<-grepl("ertoli",rawData,fixed=TRUE)

algorithm<-algName(F,F)

cData<-add(rawData,paste0("_",algorithm,"_confused"))
type<-atype(cData)

dir1<-paste0(dir)

stageGeneInteractionModel<-F
stepAIC<-F

for(mf in c('',"rank_","confused_"))
{
  for(cData in c(rawData)) #,rawData2))
  {
    combi<-cData==rawData
    sertoli<-grepl("ertoli",rawData,fixed=TRUE)
    
    algorithm<-algName(stepAIC,stageGeneInteractionModel)
    
    cData<-add(cData,paste0("_",algorithm,ifelse(combi,"_confused",ifelse(sertoli,"_sism","_sisf"))))
    type<-atype(cData)
    
    if(!exists("d")) d<-readRDS(cData)    
    
    d$stage<-stager(d$stage)
    
    dir1<-paste0(dir)
    dummy<-R.utils::mkdirs(dir1,mustWork=TRUE)
    
    for(score in c("sisf","sism","sismf"))    
    {
      if((mf=='confused_')&&(score!='sismf')) next
         
      for(genotype in c("",unique(d$genotype)))
      {
        if(genotype=='')
        {
          d0<-d
        } else
        {
          d0<-d[,d$genotype==genotype]
        }
        
        message(paste(mf,score,genotype,sep='\t'))
      
        if(mf=='')
        {
          suppressMessages(ggsave(filename = paste0(dir1,mf,score,"_",genotype,'.png'),
                plot_cells(d0,color_cells_by=score,show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                          ggtitle(paste(score,genotype))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
        } else if(mf=='rank_')
        {
          suppressMessages(ggsave(filename = paste0(dir1,mf,score,"_",genotype,'.png'),
                 plot_cells(d0,color_cells_by=paste0(score,'_rank'),show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                         ggtitle(paste(score,"rank",genotype))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
        } else
        {
          suppressMessages(ggsave(filename = paste0(dir1,mf,score,"_",genotype,'.png'),
                                  plot_cells(d0[,d0$confused],color_cells_by='sismf_rank',show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                                    ggtitle(paste(score,"rank",genotype))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
        }
        
        for(stage in unique(d0$stage))
        {
          sdir<-paste(dir1,"stage/")
          dummy<-R.utils::mkdirs(sdir,mustWork=TRUE)
          
          d1<-d0[,d0$stage==stage]
          
          if(mf=='')
          {
            suppressMessages(ggsave(filename = paste0(sdir,mf,score,"_",genotype,'_',stage,'.png'),
                   plot_cells(d1,color_cells_by=score,show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                     ggtitle(paste(score,genotype,stage))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
          } else if(mf=='rank_')
          {
            suppressMessages(ggsave(filename = paste0(sdir,mf,score,"_",genotype,'_',stage,'.png'),
                   plot_cells(d1,color_cells_by=paste0(score,'_rank'),show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                     ggtitle(paste(score,"rank",genotype,stage))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
          } else
          {
            suppressMessages(ggsave(filename = paste0(sdir,mf,score,"_",genotype,'_',stage,'.png'),
                                    plot_cells(d0[,d0$confused],color_cells_by='sismf_rank',show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                                      ggtitle(paste(score,"rank",genotype,stage))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
          }          
          
        }
      }
      
      for(label in unique(d$indiv_label))
      {
        ldir<-paste(dir1,"samples/")
        dummy<-R.utils::mkdirs(ldir,mustWork=TRUE)
        
        d1<-d0[,d0$indiv_label==label]
        
        if(mf=='')
        {
          suppressMessages(ggsave(filename = paste0(ldir,mf,score,"_",genotype,'_',stage,'.png'),
                 plot_cells(d1,color_cells_by=score,show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                   ggtitle(paste(score,genotype,stage))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
        } else if(mf=='rank_')
        {
          suppressMessages(ggsave(filename = paste0(ldir,mf,score,"_",genotype,'_',stage,'.png'),
                 plot_cells(d1,color_cells_by=paste0(score,'_rank'),show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                   ggtitle(paste(score,"rank",genotype,stage))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
        } else
        {
          suppressMessages(ggsave(filename = paste0(ldir,mf,score,"_",genotype,'_',stage,'.png'),
                                  plot_cells(d0[,d0$confused],color_cells_by='sismf_rank',show_trajectory_graph=F,label_cell_groups=T,label_branch_points=F,label_root=F,label_leaves=F)+
                                    ggtitle(paste(score,"rank",genotype,stage))+ylim(-15,15)+xlim(-15,15),dpi=dpi))
        }          
      }
    }
  }
}
