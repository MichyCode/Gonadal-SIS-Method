options(error=traceback) 

require(dplyr,quietly = TRUE,warn.conflicts=FALSE)
suppressMessages(require(pROC,quietly = TRUE,warn.conflicts=FALSE))
require('RMySQL',quietly=TRUE,warn.conflicts=FALSE)
require(progress,quietly=TRUE,warn.conflicts=FALSE)

require(doParallel,quietly=TRUE,warn.conflicts=FALSE)
registerDoParallel(4)

source("~/R/sql/sql.r")
source('znrf3.r')

grab_en_id <- function(gene, cds){
  a  <- match(gene, fData(cds)$gene_short_name)
  x <- fData(cds)$id[a]
  return(x)
}

#expressions<-c(5,10,20,50,100)

# Create algorithm name

stageGeneInteractionModel<-FALSE
stepAIC<-FALSE
algorithm<-algName(stepAIC,stageGeneInteractionModel)

message(paste0("Algorithm: ",algorithm))

#

for(rd in c(rawData)) #,rawData2)) # Loop through combines, xx and xy datasets
{
  genotypes<-getGenotypes(rd)

  #
  
  if(exists('cds')) rm(cds) # Make sure old data has gone
  
  ld<-add(rd,"_long")

  dtype<-atype(ld)

  message(paste0("\tDataset:\t",dtype,"\t",ld,"\tWith genotypes:\t",ifelse("Sertoli" %in% genotypes,"sertoli\t",""),ifelse("Granulosa" %in% genotypes,"granulosa\t","")))
  
  # Load data

  if(exists('d')) rm(d)
  
  #cromox<-c('X','Y')

	for(experiment in genotypes) # Experiment name - what method we used to create SIS identifier gene list
  {
	  type<-atype(experiment)
	  
	  if((type!='all')&&(tolower(experiment)!=tolower(type))) next

	  message(paste0("\t\tCell type\t",experiment))
	  
	  # Make a name with experiment+missing chromasomes

	  experiment2<-experiment

	  for(c in c('X','Y')) experiment2<-paste0(experiment2,'_',c);

	  targetcacheFile<-modelFile(ld,stageGeneInteractionModel,stepAIC,experiment)
	  
	  gm<-NA
	  sm<-NA
	  
	  if(file.exists(targetcacheFile))
	  {
	    message(paste0("\t\t\tFinal model found: ",targetcacheFile))
	    
	    message(paste0("\t\t\t\tLoading cellset data - ",targetcacheFile))
	    cellset<-readRDS(file=targetcacheFile)
	    
	    hasGranulosa<-"xx" %in% unique(cellset$genotype)
	    hasSertoli<-"xy" %in% unique(cellset$genotype)
	    
	    if((experiment=='Granulosa')&&!hasGranulosa) message("No granulosa")
	    if((experiment=='Sertoli')&&!hasSertoli) message("No sertoli")
	    
	    # Model file
  
	    mFile<-paste0('/scratch/deep-znrf3-finalfinal/methods/',experiment,'_X_Y_svmLinear.rds')

	    message(paste0("\t\t\t\tLoading fitted model - ",mFile));

	    m<-readRDS(file=mFile)
	    
	    if(!exists('cds'))
	    {
	      message(paste0("\t\t\t\tLoading raw data set - ",rd))    	      
	      cds<-readRDS(rd)
	    }
	    
	    # Add sis scores to monocle data object

	    genes<-names(m$coefficients[4:length(m$coefficients)])
	    norm_counts <- normalized_counts(cds,norm_method="size_only")
	    
	    if(experiment=='Sertoli')
	    {
	      sm<-m
	      
	      sis<-predict(m,newdata=cellset)
	      pData(cds)$sism<-sis*100; # pmax(0,pmin(100,sis*100))
	      mode<-'sism'
	      
	      pData(cds)$SgeneX5<-0

	      for(gene in genes)
	      {
	        gene1 <- grab_en_id(gene, cds)
	        
	        if(!is.na(gene1))
	        {
	          g1 <- norm_counts[gene1, ] # extract gene row.	        
	          cds$SgeneX5<-cds$SgeneX5+(g1>=5)
	        }
	      }
	      
	      pData(cds)$SgeneX50<-0
	      
	      for(gene in genes)
	      {
	        gene1 <- grab_en_id(gene, cds)
	        
	        if(!is.na(gene1))
	        {
	          g1 <- norm_counts[gene1, ] # extract gene row.	        
	          cds$SgeneX50<-cds$SgeneX50+(g1>=50)
	        }
	      }	      
	      
	    }	else
	    {
	      sm<-m
	      
	      sis<-predict(m,newdata=cellset)
	      pData(cds)$sisf<-sis*100; #pmax(0,pmin(100,sis*100))
	      mode<-'sisf'
	      
        #save(file='/scratch/sisf.rData',list=ls(all.names=TRUE))
	      
	      pData(cds)$GgeneX5<-0
	      
	      for(gene in genes)
	      {
	        gene1 <- grab_en_id(gene, cds)
	        
	        if(!is.na(gene1))
	        {
	          g1 <- norm_counts[gene1, ] # extract gene row.	        
	          cds$GgeneX5<-cds$GgeneX5+(g1>=5)
	        }
	      }
	      
	      pData(cds)$GgeneX50<-0
	      
	      for(gene in genes)
	      {
	        gene1 <- grab_en_id(gene, cds)
	        
	        if(!is.na(gene1))
	        {
	          g1 <- norm_counts[gene1, ] # extract gene row.	        
	          cds$GgeneX50<-cds$GgeneX50+(g1>=50)
	        }
	      }	      
	    }      	      
	    
	    if((length(cds$sisf)>0)&&(length(cds$sism)>0))
	    {
	      for(gene in genes)
	      {
	        gene1 <- grab_en_id(gene, cds)
	        
	        if(!is.na(gene1))
	        {
	          g1 <- norm_counts[gene1, ] # extract gene row.	        
	          cds$GgeneX5<-cds$GgeneX5+(g1>=5)
	        }
	      }
	      
	      #cell_types<-unique(cds$cell_type)
      
        #g<-c('Sertoli','Sertoli_Cycling','Pre-Sertoli','Pre-Supporting_1','Granulosa',"Pre-Supporting_3")
	      
	      #stop("FDR done")
	      
	      cutoff<-0.05

	      # Right hand box
	      cutoffgx<-quantile(cds[,(cds$genotype=='xx')&(cds$stage=='E14.5')&(cds$cell_type=='Granulosa')]$sism,1-cutoff/2)
	      
	      # Left hand box
	      
	      cutoffsx<-quantile(cds[,(cds$genotype=='xy')&(cds$stage=='E14.5')&(cds$cell_type=='Sertoli')]$sism,cutoff/2)
	      cutoffsy<-quantile(cds[,(cds$genotype=='xy')&(cds$stage=='E14.5')&(cds$cell_type=='Sertoli')]$sisf,1-cutoff/2)
	      
	      message(paste("Cutoffs:\t",cutoffgx,cutoffsx,cutoffsy));
	      
	      #pData(cds)$confused0<-(cds$sism>=cutoffs)&(cds$sisf>=cutoffg)
	      
  	    pData(cds)$sismf<-cds$sism*cds$sisf
  	    
  	    pData(cds)$sism_rank<-perc.rank2(cds[,cds$genotype!='xy_hom']$sism,cds$sism)
  	    pData(cds)$sisf_rank<-perc.rank2(cds[,cds$genotype!='xy_hom']$sisf,cds$sisf)
	      pData(cds)$sismf_rank<-perc.rank2(cds[,cds$genotype!='xy_hom']$sism*cds[,cds$genotype!='xy_hom']$sisf,cds$sisf*cds$sism)
	      
	      pData(cds)$confused_orig<-abs(cds$sism-cds$sisf)<sqrt(cds$sismf)
	      pData(cds)$confused<-(cds$sism>=cutoffgx)&((cds$sism<=cutoffsx)|(cds$sisf>=cutoffsy))
	      
	      pData(cds)$coex<-abs(cds$sism-cds$sisf)
	      
	      #pData(cds)$confused1<-(cds$sism>=0.5)&(cds$sisf>=0.5)
	      
	      
	      #cutm<-quantile(cds[,cds$genotype=='xx']$sism,c(cut))
	      #cutf<-quantile(cds[,cds$genotype=='xy']$sisf,c(cut))
	      
	      #pData(cds)$confused<-(cds$sism>=cutm)&(cds$sisf>=cutf)
	      
	      mode<-'confused'
	      
	      #save(file='/scratch/sismf.rData',list=ls(all.names=TRUE))    	      
	    }
	    
	    message(paste(dtype,mode))
	    
	    if(((dtype=='all')&&(mode=='confused'))||((dtype=='Granulosa')&&(mode=='sisf'))||((dtype=='Sertoli')&&(mode=='sism')))
	    {
	      newFile<-add(rd,paste0("_",'svm',"_",mode))
	    
	      message(paste0(targetcacheFile,"\n",newFile,"\n",older(targetcacheFile,newFile))) 
	      
	      if(older(targetcacheFile,newFile)) # Only save if file missing, or model file newer than existing save
	      {
	        message(paste0("\t\t\t\tSaving monocle object with sis scores - ",newFile))
	        
	        saveRDS(file=newFile,object=cds)
	      }
	    }
	  } else
	  {
	    message(paste0("Error - cache file not found - ",targetcacheFile))
	  }
  }
}
