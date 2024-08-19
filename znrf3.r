###############################################################################
# Version is the primary key, always increment. Project name can be anything. #
###############################################################################

version<-'finalfinal'              # <- Change this
project<-'deep-znrf3'
db<-'single_cell'
resultsDir<-'/scratch/'

#interesting<-c('Sertoli','Granulosa','Pre-Supporting_1','Pre-Supporting_2','Pre-Supporting_3','Pre-Sertoli') # 'Pre-Supporting_Bmyc',
#interesting<-c('Cycling_Supporting','Pre_Supporting_2','Pre_Supporting_1','Granulosa','SLC','Pre_Sertoli','Sertoli')
interesting<-c('Cycling-Supporting','Pre-Supporting_2','Pre-Supporting_1','Granulosa','SLC','Pre-Sertoli','Sertoli')

message(paste0("****** ",project,"-",version," ",db," ",resultsDir," ******"))

# Graph stuff

require(ggplot2)

source('../zcolours.r')

###############################################################################

options(error=traceback)

require(dplyr,quietly = TRUE,warn.conflicts=FALSE)
suppressMessages(require(monocle3,quietly = TRUE,warn.conflicts=FALSE))

mtime<-function(f)
{
  as.numeric(file.info(f)$mtime)
}

perc.rank <- function(x) trunc(rank(x))/length(x)
perc.rank2 <- function(x, xo)
{
  zx<-c()
  
  for(z in xo)
  {
    p<-length(x[x <= z])/length(x)
    zx<-c(zx,p)
    
    #message(paste(length(zx),p,sep='\t'))
  }
  
  zx
}

perc.rank3 <- function(x, xo)
{
  s<-length(x)
  
  length(x[x <= xo])/s
}

algName<-function(stepAIC,stageGeneInteractionModel)
{
  algorithm<-'';
  
  if(stepAIC) algorithm<-paste0(algorithm,"Aic")
  
  if(stageGeneInteractionModel) algorithm<-paste0(algorithm,ifelse(stringr::str_length(algorithm)>0,"_",""),"SGIM")
  
  algorithm<-ifelse(stringr::str_length(algorithm)==0,"Normal",algorithm)
  
  algorithm
}

getGenotypes<-function(f)
{
  genotypes<-c()
  
  if(grepl("sertoli",rd, fixed = TRUE)||grepl("Sertoli",rd, fixed = TRUE)||grepl("xy",rd, fixed = TRUE)) genotypes<-c(genotypes,"Sertoli")
  if(grepl("granulosa",rd, fixed = TRUE)||grepl("Granulosa",rd, fixed = TRUE)||grepl("xx",rd, fixed = TRUE)) genotypes<-c(genotypes,"Granulosa")
  
  if(length(genotypes)==0) genotypes<-c("Sertoli","Granulosa")
  
  genotypes;
}

stager<-function(stage)
{
  paste0(substring(stage,2)," dpc")
}

sstop<-function(m="",rc=0)
{
  if(interactive())
  {
    stop(m)
  } else
  {
    message(m)
    quit(save="no",status = rc)  
  }
}

older<-function(filea,fileb)
{
  if(!file.exists(fileb))
  {
    T
  } else
  {
    agea<-mtime(filea)
    ageb<-mtime(fileb)
  
    #message(paste(filea,agea,fileb,ageb))
    
    ageb<agea
  }
}

add<-function(s,extra) 
{
  f<-tools::file_path_sans_ext(s)
  e<-tools::file_ext(s)
  
  paste0(f,extra,'.',e)
}

cellify<-function(cellnames,samplenames)
{
  stringr::str_replace(cellnames,"-.*",paste0("_",samplenames))
}

atype<-function(ld)
{
  if(grepl("sertoli",ld, fixed = TRUE)||grepl("Sertoli",ld, fixed = TRUE))
  {
    type<-'Sertoli';
  } else if(grepl("granulosa",ld, fixed = TRUE)||grepl("Granulosa",ld, fixed = TRUE))
  {
    type<-'Granulosa';
  } else
  {
    type<-'all'
  }
  
  type
}

`%notin%` <- Negate(`%in%`)

# Identify saved model filename

modelFile<-function(ld,stageGeneInteractionModel,stepAIC,experiment)
{
  dtype<-atype(ld)
  
  algorithm<-'';
  
  if(stepAIC) algorithm<-paste0(algorithm,"Aic")
  
  if(stageGeneInteractionModel) algorithm<-paste0(algorithm,ifelse(stringr::str_length(algorithm)>0,"_",""),"SGIM")
  
  algorithm<-ifelse(stringr::str_length(algorithm)==0,"Normal",algorithm)
  
  type<-atype(experiment)
  
  experiment2<-experiment
  
  for(c in c('X','Y')) experiment2<-paste0(experiment2,'_',c);
  
  paste0(modelDir,'sex_gene_',dtype,"_",algorithm,"_all_",experiment2,'.rds')
}

ssaveRDS<-function(file,object)
{
  saveRDS(file=file,object=object,compress=F)
  system(paste("pxz -9 ",file))
  file.rename(from=paste0(file,".xz"),to=file)
}

#####################
# Project constants #
#####################

# Database tables

geneDB<-paste0(db,'.gene_roc',version)
modelDB<-paste0(db,'.model',version)
methodDB<-paste0(db,'.method',version)

# Results data

dir<-paste0(resultsDir,project,"-",version)

dataDir<-paste0(dir,"/data/")
cellsetDir<-paste0(dir,"/cellset/")
modelDir<-paste0(dir,"/models/")
methodDir<-paste0(dir,"/methods/")
graph<-paste0(dir,"/graphs/")

if(version=='8')
{
  rawData<-paste0(dataDir,"znrf3_exp9_all_obj_cell_type_k17_labelled.rds")
  rawData2<-c(paste0(dataDir,"znrf3_xx_wt_labelled_obj.rds"),paste0(dataDir,"znrf3_xy_wt_labelled_obj.rds"))
} else if(version=='8a')
{
  rawData<-paste0(dataDir,"znrf3_exp9_all_obj_cell_type_k15_labelled.rds")
  rawData2<-c(paste0(dataDir,"znrf3_xx_wt_labelled_obj.rds"),paste0(dataDir,"znrf3_xy_wt_labelled_obj.rds"))
} else if (version=='7a')
{
  rawData<-paste0(dataDir,"znrf3_obj_all_gens_clusters_labelled.rds")
  rawData2<-c(paste0(dataDir,"znrf3_xx_wt_all_tps_granulosa_lab_7-10-2020.rds"),paste0(dataDir,"znrf3_xy_wt_all_tps_sertoli_lab_7-10-2020.rds"))
} else if(version=='final')
{
  rawData<-paste0(dataDir,"znrf3_all_full_obj_labelled.rds")
  rawData2<-c(paste0(dataDir,"znrf3_xx_wt_clustered_k13.rds"),paste0(dataDir,"znrf3_xy_wt_clustered_k7.rds"))
} else if(version=='final4')
{
  rawData<-paste0(dataDir,"znrf3_obj_all_gens_clusters_labelled.rds")
  rawData2<-c(paste0(dataDir,"xx.rds"),paste0(dataDir,"xy.rds"))
} else if(version=='finalfinal')
{
  rawData<-paste0(dataDir,"all_wt_obj_pc55_k19.rds")
  rawData2<-c(paste0(dataDir,"xx_wt_obj_pc25_k8.rds"),paste0(dataDir,"znrf3_xy_wt_clustered_k7.rds"))
}

longData<-paste0(dataDir,"_long.rds")

cData<-paste0(dataDir,"dummy_with_confused.rds")
cData2<-add(rawData2,"with_sis")

gdir<-paste0(dir,'/graphs/')

houseKeepingData<-paste0(dataDir,'Housekeeping_Genes_Mouse.RData')

dummy<-R.utils::mkdirs(dataDir,mustWork=T)
dummy<-R.utils::mkdirs(cellsetDir,mustWork=T)
dummy<-R.utils::mkdirs(modelDir,mustWork=T)
dummy<-R.utils::mkdirs(methodDir,mustWork=T)
dummy<-R.utils::mkdirs(graph,mustWork=T)

# Gonad development genes, from Riassa:

ovaryDevGenes<-c('Wnt4','Fst','Foxl2','Ctnnb1','Rspo1')
testesDevGenes<-c('Sry','Sox9','Amh','Fgf9','Pgd2','Wt1','Map3k4','Gata4','Igf1r','Cbx2','Dmty1','Nr5a1')

# for (obj in ls()) { message(obj); print(object.size(get(obj)), units='auto') }

#options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

# get Ensemble ID
grab_en_id <- function(gene, cds){
  a  <- match(gene, fData(cds)$gene_short_name)
  x <- fData(cds)$id[a]
  return(x)
}

require(ggpubr)

youLegend <- function(plot,file,width=NA,height=NA)
{
	ggsave(filename=file,plot+theme(legend.position ="none"),dpi=dpi,width=width,height=height)

	tryCatch({
		lfile<-add(file,"_legend")

		leg<-get_legend(plot)
		ggsave(filename=lfile,as_ggplot(leg),dpi=dpi)
	},error=function(e)
      {
      })
}
