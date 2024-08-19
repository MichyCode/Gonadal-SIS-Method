require(dplyr,quietly = TRUE,warn.conflicts=FALSE)
require(pROC,quietly = TRUE,warn.conflicts=FALSE)
require('RMySQL',quietly=TRUE,warn.conflicts=FALSE)
require(progress,quietly=TRUE,warn.conflicts=FALSE)

require(doParallel,quietly=TRUE,warn.conflicts=FALSE)
registerDoParallel(4)

source("~/R/sql/sql.r")
source('znrf3.r')

# To re-run delete files in 'models' directory

for(stageGeneInteractionModel in c(FALSE))
{
  for(stepAIC in c(FALSE))
  {
    algorithm<-algName(stepAIC,stageGeneInteractionModel)
    
    message(paste0("Algorithm: ",algorithm))

    for(rd in c(rawData,rawData2)) # Loop through combines, xx and xy datasets
    {
      genotypes<-getGenotypes(rd)
      
      if(exists('cds')) rm(cds)
      
      ld<-add(rd,"_long")
    
      type<-atype(ld)
    
      message(paste0("\tDataset:\t",type,"\t",ld))
      
      # Update q values for genes
    
      if(singleQuery("SELECT COUNT(*) FROM ",geneDB," WHERE (q_x IS NULL or q_intercept IS NULL) AND dataset='",ld,"'")>0)
      {
        message(paste0("**** Updating gene q values for '",type,"'"))
      
        ps<-tableQuery("SELECT bod_id,p_x,p_intercept FROM ",geneDB," WHERE dataset='",ld,"'")
    
        qx<-qvalue::qvalue(ps$p_x)$qvalues
        qi<-qvalue::qvalue(ps$p_intercept)$qvalues
    
        for(i in 1:nrow(ps))
        {
          executeQuery("UPDATE ",geneDB," SET q_x=",qx[i],",q_intercept=",qi[i]," WHERE bod_id=",ps$bod_id[i]," AND  dataset='",ld,"'")
        }
      }
    
      R.utils::mkdirs(modelDir,mustWork=TRUE)
    
      # Load data
    
      if(exists('d')) rm(d)
      
      cromox<-c('X','Y')
    
    	for(experiment in genotypes) # Experiment name - what method we used to create SIS identifier gene list
      #foreach(experiment=genotypes) %dopar%
      {
    	  type<-atype(experiment)
    	  
    	  if((type!='all')&&(experiment!=type)) next
    
    	  message(paste0("\t\tCell type\t",experiment))
    	  
    	  # Make a name with experiment+missing chromasomes
    
    	  experiment2<-experiment
    
    	  for(c in c('X','Y')) experiment2<-paste0(experiment2,'_',c);
    
    	  targetcacheFile<-modelFile(ld,stageGeneInteractionModel,stepAIC,experiment)
    	  
    	  if(!older(ld,targetcacheFile)&&FALSE)
    	  {
    	     message(paste0("\t\t\tFinal model found: ",targetcacheFile))
    	  } 
    	  else
    	  {
    	    message(paste0("\t\t\tBuilding model ",ld,"\t->\t",targetcacheFile,ifelse(file.exists(targetcacheFile),"\t(Missing)","\t(out of date)")))
    	    
      	  # Load candidate genes - if present
      	  
      	  #bsql<-paste0("SELECT * FROM ",geneDB," WHERE model='sis~x+stage' AND q_x<0.05 and q_intercept<0.05 AND exp='",experiment,"' AND roc>=0.7 AND dataset='",ld,"' ORDER BY roc DESC")
    	    bsql<-paste0("SELECT * FROM ",geneDB," WHERE model='sis~x+stage' AND exp='",experiment,"' AND roc>=0.7 AND dataset='",ld,"' ORDER BY roc DESC")
      	  best<-tableQuery(bsql) # Get best SIS identifier genes
      	  
      	  if(nrow(best)==0) # No candidates?
      	  {
      	    message(paste0("\t\t\t**** No candidate genes found"))
      	  } else
      	  {
        	  # Make dataset
        	  
        	  if(!exists("d"))
        	  {
        	    message(paste0("\t\t\tLoading long data from '",ld,"'"))
        	    d<-readRDS(file=ld)
        	  }
        	  
        		cellset<-distinct(d,d$cell,.keep_all=TRUE) # Take 1 of each cell (doesn't matter which)
        	  #cellset<-cellset[,-c(1,2,11,12,13,14,15)]            # Remove added stuff from 'distinct' function
        
        	  # Set up the training set to be the correct wt, and the training SIS to be the correct cell typ
        	  
        		if(experiment=='Granulosa')
        	  {
        		  cellset$sis<-as.integer(cellset$cell_type_g) # Create a SIS score
        	    training<-cellset$training_g
        	  } else if(experiment=='Sertoli')
        	  {
        		  cellset$sis<-as.integer(cellset$cell_type_s) # Create a SIS score
        	    training<-cellset$training_s
        	  }
        
        	  #message(paste("\tBuilding models for ",experiment2))
    
        	  
        	  model<-'' # Initial model with no genes
        
        	  index<-1 # Number of genes in model
        
        	  for(i in 1:nrow(best)) # Loop through the genes (need index for saving cache files)
        	  {
        	    #message(paste("\t\tTrying gene ",i))
        	    
          		gene<-best[i,'gene'];
        		  chromosome<-singleQuery("SELECT chromosome FROM mgi_cache.marker WHERE symbol='",gene,"'");
        
        		  if(chromosome %in% c('X','Y')) next;
        
      			  # Data in d0 is sparse so we need to check it exists and act accordingly
      			  
         	 	  d0<-d[d$gene==gene,]
      
         	 		x2<-d0[d0$cell %in% cellset[,'cell'],'x'] # Get available data
         	 		x<-c()
         	 		
         	 		for(z in (cellset[,'cell']  %in% d0$cell))
         	 		{
         	 		  if(z) # data exists, so add from x2
         	 		  {
         	 		    x<-c(x,x2[1])
         	 		    x2<-x2[-1]
         	 		  }
         	 		  else { # No data, so assume 0
         	 		    x<-c(x,0)
         	 		  }
         	 		}
         	 		
         	 		#message("\t\t\tAdding new gene data to cellset")
         	 		
        		  cellset<-cbind(cellset,x)
          
        		  if(grepl("^[[:digit:]]+",gene)) gene<-paste0("gene_",gene)
      
          	  gene<-stringr::str_replace(gene,"-","_")
          
        		  colnames(cellset)[ncol(cellset)]<-gene
    
        			#message("\t\t\t\tRunning model...");
        			
        		  if(grepl("^[[:digit:]]+",gene)) gene<-paste0("gene_",gene)
          
        		  gene<-stringr::str_replace(gene,"-","_")
          
        		  oldmodel<-model
        		  
        		  stages<-length(unique(cellset$stage))>1
        		  
      		    model<-paste0(model,ifelse(stringr::str_length(model)==0,"","+"),gene)

   		        formula<-paste0("sis~",ifelse(stages,"stage","1"),"+",model)

        	    m<-stats::glm(as.formula(formula),subset=training,data=cellset) # Fit a model
      		    s<-summary(m)                                      # Get the summary

        		  m$data<-NULL
        		  
        		  cellset0<-cellset
          	  cellset0$p<-predict(m,newdata = cellset0);
          	  
        			cellset0<-cellset0[training,];
  
        			r<-roc(cellset0$sis>0.5,cellset0$p,quiet=TRUE); # ROC curve
        		  a<-auc(r);                                      # And the ROC curve AUC
          	  
          	  #
          	  
      	      message(paste("\t\t\t",i,gene,best[i,'roc'],a,index));
      	        
      			  executeQuery("REPLACE INTO ",modelDB," (alg,cell_type,gene,genes,added_roc,model_roc,dataset) VALUES (",
    	  		                "'",algorithm,"','",experiment2,"','",gene,"',",index,",",best[i,'roc'],",",a,",'",ld,"')")
    
        	    if(index==50)
        	    {
        	      message(paste0("\t\tSaving: ",targetcacheFile))
    ssaveRDS(file=targetcacheFile,object=cellset) # Save model
        	      message(paste0("\t\tSaving: ",add(targetcacheFile,'_model')))
    ssaveRDS(file=add(targetcacheFile,'_model'),object=m) # Save model
        	      
        	      break;
        	    }
        	    
        	    index<-index+1
        	  }
      	  }
    	  }
    	}
    }
  }
}
