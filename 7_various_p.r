options(error=traceback)

rm(list=ls())

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

if(!file.exists(cData)) sstop(paste(cData,"does not exist"),1)

if(!exists("d"))
{
  message(paste0("\tLoading data from ",cData))

  d<-readRDS(cData)

  z<-data.frame(d$sism,d$sisf,d$sismf,d$stage,d$cell_type,d$genotype,d$label,d$confused,colnames(d))
  colnames(z)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","sample","confused","cell")

  z$stage<-stager(as.factor(z$stage))
  z$genotype<-as.factor(z$genotype)
  
  #z[is.na(z$confused_orig),'confused_orig']<-FALSE

  data<-readRDS("/scratch/deep-cd1-final/data/cd1b6j_monocle3_newcelldataset_Normal_confused.rds");

  zcd1<-data.frame(data$sism,data$sisf,data$sismf,data$time_ok_ok,data$predicted_annot_B6,data$genotype_ok,data$indiv_label,data$confused,colnames(data))
  colnames(zcd1)<-c("sis_m","sis_f","sis_mf","stage","cell_type","genotype","confused","cell")
   
  zcd1$stage<-stager(as.factor(zcd1$stage))
  zcd1$genotype<-as.factor(tolower(zcd1$genotype))
  zcd1genotype<-factor(zcd1$genotype,level=c('xx','xy'))
}

z1<-z[z$cell_type %in% interesting,]
z1$group<-as.factor(ifelse(z1$cell_type %in% c('Pre-Sertoli','Sertoli'),"Sgroup","Ggroup"))

# Stage

z1<-z1[z1$stage=='12.5 dpc',]
zcd1<-zcd1[zcd1$stage=='12.5 dpc',]

for(genotype in c('xy'))
{
	message("")
	message(genotype)
	message("")

	s<-z1[(z1$group=='Sgroup')&(z1$genotype==genotype),'sis_f']
	g<-z1[(z1$group=='Ggroup')&(z1$genotype==genotype),'sis_m']

	#

	zcd11<-zcd1[zcd1$cell_type %in% interesting,]
	zcd11$group<-as.factor(ifelse(zcd11$cell_type %in% c('Pre-Sertoli','Sertoli'),"Sgroup","Ggroup"))

	sc<-zcd11[(zcd11$group=='Sgroup')&(zcd11$genotype==genotype),'sis_f']
	gc<-zcd11[(zcd11$group=='Ggroup')&(zcd11$genotype==genotype),'sis_m']

	cat("S group:\t",paste(unique(z1[z1$group=='Sgroup','cell_type']),sep='\t'),"\n")
	cat("G group:\t",paste(unique(z1[z1$group=='Ggroup','cell_type']),sep='\t'),"\n")

	cat("S cd1 group:\t",paste(unique(zcd11[zcd11$group=='Sgroup','cell_type']),sep='\t'),"\n")
	cat("G cd1 group:\t",paste(unique(zcd11[zcd11$group=='Ggroup','cell_type']),sep='\t'),"\n")


	message("Mean(s group)=",mean(s),"\tvar(s group)=",var(s),"\tn=",length(s))
	message("Mean(g group)=",mean(g),"\tvar(g group)=",var(g),"\tn=",length(g))
	message("Mean(s cd1 group)=",mean(sc),"\tvar(s cd1 group)=",var(sc),"\tn=",length(sc))
	message("Mean(g cd1 group)=",mean(gc),"\tvar(g cd1 group)=",var(gc),"\tn=",length(gc))

	message("")

        if((length(s)>0)&&(length(g)>0))
	{
		message("B6N g v s:");
		message("p(mean)=",wilcox.test(s,g)$p.value)
		message("p(var)=",var.test(s,g)$p.value)
	}

        if((length(sc)>0)&&(length(gc)>0))
	{
		message("CD1 g v s:");
		message("p(mean)=",wilcox.test(sc,gc)$p.value)
		message("p(var)=",var.test(sc,gc)$p.value)
	}

        if((length(sc)>0)&&(length(s)>0))
	{
		message("S CD1 v B6N");
		message("p(mean)=",wilcox.test(sc,s)$p.value)
		message("p(var)=",var.test(sc,s)$p.value)
	}

        if((length(gc)>0)&&(length(g)>0))
	{
		message("G CD1 v B6N");
		message("p(mean)=",wilcox.test(gc,g)$p.value)
		message("p(var)=",var.test(gc,g)$p.value)
	}
}
