require(monocle3,quietly=TRUE,warn.conflicts=FALSE)
require(dplyr,quietly=TRUE,warn.conflicts=FALSE)
require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)
require(ggthemes,quietly=TRUE,warn.conflicts=FALSE)
require(gridExtra,quietly=TRUE,warn.conflicts=FALSE)
require(reshape2,quietly=TRUE,warn.conflicts=FALSE)

source('znrf3.r')

stageGeneInteractionModel<-F
stepAIC<-F
cData<-rawData

combi<-cData==rawData
sertoli<-grepl("ertoli",rawData,fixed=TRUE)||grepl("xy_wt",rawData,fixed=TRUE)

algorithm<-algName(stepAIC,stageGeneInteractionModel)

cData<-add(cData,paste0("_",algorithm,ifelse(combi,"_confused",ifelse(sertoli,"_sism","_sisf"))))
type<-atype(cData)

graph<-paste0(graph,"cutoffs")
dummy<-R.utils::mkdirs(graph,mustWork=T)

message(paste0("Algorithm: ",algorithm,'\t','Dataset type: ',type))

gdir0<-paste0(gdir,algorithm,"/",type,"/")

dummy<-R.utils::mkdirs(gdir0,mustWork=TRUE)

cData<-'/scratch/deep-znrf3-finalfinal/data/all_wt_obj_pc55_k19_svm_confused.rds'

if(!file.exists(cData)) sstop(paste(cData,"does not exist"),1)

if(!exists("d"))
{
  message(paste0("\tLoading data from ",cData))

  d<-readRDS(cData)

  z<-data.frame(d$sism,d$sisf,d$sismf,d$stage,d$cell_type,d$genotype,d$label,d$confused,d$confused_orig,colnames(d))
  colnames(z)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","sample","confused","confused_orig","cell")

  z$stage<-as.factor(z$stage)
  z$genotype<-as.factor(z$genotype)
  
  z[is.na(z$confused_orig),'confused_orig']<-FALSE
}

# Stats

interesting<-c('Sertoli','Granulosa','Pre-Supporting 1','Pre-Supporting 2','Pre-Supporting 3','Pre-Supporting Bmyc','Pre-Sertoli')

z1<-z[z$stage=='E12.5',];

#graph<-paste0(graph,"/components/")

pl<-list()

for(cutmode in c('global_cut')) #,'sg_cut'))
{
  for(cut in c(.95)) # ,.975,.99)) # Graph cutoff v numbers v P
  {
    if(cutmode=='global_cut')
    {
      cutm<-quantile(z1[z1$genotype=='xx','sis_m'],c(cut))
      cutf<-quantile(z1[z1$genotype=='xy','sis_f'],c(cut))
    } else if (cutmode=='sg_cut')
    {
      cutm<-quantile(z1[(z1$genotype=='xx')&(z1$cell_type=='Granulosa'),'sis_m'],c(cut))
      cutf<-quantile(z1[(z1$genotype=='xy')&(z1$cell_type=='Sertoli'),'sis_f'],c(cut))
    }
  
    message(cutm)
    message(cutf)
    #stop("!")
    
    qs<-NA
    
    celltypes<-unique(z1$cell_type);
    
    for(pass in 1:3)
    {
      ps<-c()
    
      for(ci in 1:length(celltypes))
      {
        celltype<-celltypes[ci]
        
        z2<-z1[z1$cell_type==celltype,] 
        
        z2$con<-(z2$sis_m>=cutm)&(z2$sis_f>=cutf)
        #z2$con<-(z2$sis_m>=50)&(z2$sis_f>=50)
        
        #
        
        ldp<-sum(z2[z2$genotype=='xy_hom','confused_orig'])
        lt<-sum(!z2[z2$genotype=='xy_hom','confused_orig'])
        
        adp<-sum(z2[z2$genotype=='xy','confused_orig'])
        at<-sum(!z2[z2$genotype=='xy','confused_orig'])
        
        #
        
        t<-matrix(c(ldp,lt,adp,at),ncol=2)
        rownames(t)<-c('confused','not confused')
        colnames(t)<-c('hom','het')
        
        c<-chisq.test(t)
        
        p<-c$p.value
        if(is.na(p)) p<-1
        ps<-c(ps,p)  
        
        if(pass==2)
        {
          message(paste(cut,celltype,p,qs[ci],sep='\t'))
          
          if(celltype %in% c("Sertoli","Granulosa","Pre-Sertoli"))
          {
            for(genotype in c('xy','xx','xy_hom'))
            {
              z3<-z2[z2$genotype==genotype,]
            
              if(nrow(z3)>1)
              {
                lm_fit <- lm(sis_f ~ sis_m, data=z3)
                z3$fit<-predict(lm_fit)
              
                z4<-z3[z3$con,]
              
                if(nrow(z4)>0)
                {
                  lm_fit2 <- lm(sis_f ~ sis_m, data=z4)
                  z4$fit<-predict(lm_fit2)
                }
                
                #low<sum(z3$sis_m<)
                
                p<-ggplot(z3,aes(sis_m,sis_f)) +
                  geom_point(aes(color=con))+
                  geom_line(color='grey',data = z3,aes(x=sis_m, y=fit))+
                  geom_line(color='black',data = z4,aes(x=sis_m, y=fit))+
                  coord_cartesian(xlim=c(0,100),ylim=c(0,100))+
                
                  ggtitle(paste(celltype,"E12.5",genotype,cutmode,cut))
                
                suppressMessages(ggsave(file=paste0(graph,"/svm_ct_","_",celltype,"_",genotype,"_",cutmode,"_",cut,".png"),p))
                
              }
            }
          } 
        } else if(pass==3)
        {
          # if(qs[ci]<=0.05)
          # {
          #   message(paste(celltype))
          #   
          #   print(t)
          # }
        }
      }
      
      for(genotype in c('xy','xx','xy_hom'))
      {
        z2<-z1[z1$cell_type %in% interesting,] 
        z2$con<-(z2$sis_m>=cutm)&(z2$sis_f>=cutf)
        
        #
        
        ldp<-sum(z2[z2$genotype=='xy_hom','con'])
        lt<-sum(!z2[z2$genotype=='xy_hom','con'])
        
        adp<-sum(z2[z2$genotype=='xy','con'])
        at<-sum(!z2[z2$genotype=='xy','con'])
        
        #
        
        t<-matrix(c(ldp,lt,adp,at),ncol=2)
        rownames(t)<-c('confused','not confused')
        colnames(t)<-c('hom','het')
        
        c<-chisq.test(t)
        
        p<-c$p.value
        
        #
        
        z3<-z2[z2$genotype==genotype,]
        
        if(nrow(z3)>1)
        {
          lm_fit <- lm(sis_f ~ sis_m, data=z3)
          z3$fit<-predict(lm_fit)
          
          z4<-z3[z3$con,]
          
          if(nrow(z4)>0)
          {
            lm_fit2 <- lm(sis_f ~ sis_m, data=z4)
            z4$fit<-predict(lm_fit2)
          }
  
          p0<-ggplot(z3,aes(sis_m,sis_f)) +
            geom_point(aes(color=con))+
            geom_line(color='grey',data = z3,aes(x=sis_m, y=fit))+
            geom_line(color='black',data = z4,aes(x=sis_m, y=fit))+
            coord_cartesian(xlim=c(0,100),ylim=c(0,100))+ggtitle(paste("E12.5",genotype,cutmode,cut,"p(confused)=",p))
          
          if(genotype=='xy_hom') message(paste(genotype,cut,p))
          
          p<-ggplot(z3,aes(sis_m,sis_f)) +
            geom_point(aes(shape=con,color=cell_type))+
            geom_line(color='grey',data = z3,aes(x=sis_m, y=fit))+
            geom_line(color='black',data = z4,aes(x=sis_m, y=fit))+
            coord_cartesian(xlim=c(0,100),ylim=c(0,100))+ggtitle(paste("E12.5",genotype,cut,ifelse(genotype=='xy_hom',paste0("p(con)=",p),"")))
          
          
          
          ggsave(file=paste0(graph,"/svm_ct_","_interesting_",genotype,"_",cut,".png"),p,width = 9,height=7)
  
          if(genotype=='xx')
          {
            g1<-p0
            g4<-p
          } else if(genotype=='xy_hom')
          {
            g2<-p0
            g5<-p
          } else
          {
            g3<-p0
            g6<-p
          }
        }    
      }
      
      if(pass==1) qs<-qvalue::qvalue(ps)$qvalues
    }
    
    ggsave(file=paste0(graph,"/svm_ct_","_interesting_",cutmode,"_",cut,".png"),grid.arrange(g1,g2,g3,g4,g5,g6,nrow=2,ncol=3),width = 24,height=12)
    
    pl[[length(pl)+1]]<-g1+ theme(legend.position ="none")
    pl[[length(pl)+1]]<-g2+ theme(legend.position ="none")
    pl[[length(pl)+1]]<-g3+ theme(legend.position ="none")
    pl[[length(pl)+1]]<-g4+ theme(legend.position ="none")
    pl[[length(pl)+1]]<-g5+ theme(legend.position ="none")
    pl[[length(pl)+1]]<-g6
  }
  
  ggsave(file=paste0(graph,"/svm_ct_","_interesting_",cutmode,".png"),do.call("grid.arrange",c(pl,ncol=length(pl)/3)),width = 1920/1080*12,height=12)
}




# m<-z1[z1$genotype=='xx','sis_m']
# f<-z1[z1$genotype=='xy','sis_f']
# 
# for(dist in c("beta", "cauchy", "chi-squared", "exponential", "gamma", "geometric", "log-normal", "lognormal", "logistic", "negative binomial", "normal", "Poisson", "t""weibull"))
# {
#   message(dist)
#   #fit<-fitdistrplus::fitdistr(m,dist)
#   fit<-MASS::fitdistr(m,dist)
#   
#   message(paste(dist,fit$aic,fit$bic,sep='\t'))
# }
