#######################################################################################################################
#
# best_sex_genes
#
# This script takes the monocle data, transformed into a long format by 'extract_data', and runs simple GLM models on
# each gene to determine it's ability to classify sertoli/granulosa cells from other cell types
#
# The gene's scores are stored in a database 'single_cell' table gene_rocX
#
#######################################################################################################################

options(error=traceback)

source('znrf3.r')

require(dplyr,quietly = TRUE,warn.conflicts=FALSE)
require(pROC,quietly = TRUE,warn.conflicts=FALSE)
require('RMySQL',quietly=TRUE,warn.conflicts=FALSE)

source("~/R/sql/sql.r")

require(doParallel,quietly=TRUE,warn.conflicts=FALSE)

# To repeat this stage running the following sql or 'touch' the long data files in the data directory

# executeQuery("TRUNCATE TABLE ",geneDB);

for(rd in c(rawData)) #,rawData2[2])) # Loop through combined, xx and xy datasets
{
  ld<-add(rd,"_long")
  
  if(exists("d")) rm(d)
  
  genotypes<-getGenotypes(rd)
  
  mt<-mtime(ld) # Timestamp of data file
  
  for(experiment in genotypes) # Build dataset for Granulosa and Sertoli 
  {
    # Remove data from older version of the long data file
    
    #executeQuery("DELETE FROM ",geneDB," WHERE exp='",experiment,"' AND dataset='",ld,"' AND mtime<",mt)

    if(experiment=='Granulosa')
    {
      if(!exists("d"))      {
        message(paste0("Loading long data '",ld,"'"))
        d<-readRDS(file=ld) # Read the data
      }
      
      if(sum((d$genotype=='xx')&(d$cell_type_g))==0) next; # If no data next experiment

      d0<-d[(d$genotype=='xx'),] # Only wt
  	  d0$sis<-as.integer(d0$cell_type_g) # Create a SIS score
    } else # Sertoli
    {
      if(!exists("d"))      {
        message(paste0("Loading long data '",ld,"'"))
        d<-readRDS(file=ld) # Read the data
      }
      
      if(sum((d$genotype=='xy')&(d$cell_type_s))==0) next; # If no data next experiment
      
      d0<-d[(d$genotype=='xy'),] # Only wt
    	d0$sis<-as.integer(d0$cell_type_s) # Create a SIS score
    }

    genes<-unique(d0$gene)
    cells<-unique(d0$cell)
  
    done<-listQuery("SELECT gene FROM ",geneDB," WHERE model='sis~x+stage' AND exp='",experiment,"' AND dataset='",ld,"'");
    genes<-setdiff(genes,done)
  
    cellset<-distinct(d0,d0$cell,.keep_all=TRUE) # Take 1 of each cell (doesn't matter which)
    cellset['x']<-0;                             # Set the signal for the missing to zero
    cellset<-cellset[,-ncol(cellset)];           # Remove added stuff from 'distinct' function
  
    rm(done)

    gc(full=TRUE)
    
    memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))  
    
    message(paste0("Using ",min(as.integer(memfree/5000000),parallel::detectCores())," cores"))
    
    registerDoParallel(20); #cores=min(as.integer(memfree/5000000),parallel::detectCores()))
  
    # Loop through each gene missing from database
    
    foreach(gene=genes) %dopar%
    {
      tryCatch({
      
        d1<-d0[d0$gene==gene,]
    
        stages<-length(unique(d1$stage))>1
        
        # Add missing cells
    
        example<-cellset[!(cellset$cell %in% unique(d1$cell)),]
        
        if(nrow(example)>0)
        {
          example['gene']<-gene;               # All the same gene
          example['x']<-0;
          example['house_keeping']<-d1[1,'house_keeping']
        
          d1<-rbind(d1,example)                # Add the missing cells to the cells with the gene
        }
  
        # Analyse

        m<-stats::glm(as.formula(paste0('sis~x',ifelse(stages,'+stage',''))),data=d1) # Fit a model
        s<-summary(m)                             # Get the summary
  
        r<-roc(d1$sis>0.5,predict(m),quiet=TRUE)  # ROC curve
        a<-auc(r)                                 # And the ROC curve AUC
        effect<-m$coefficients['x']

        cat(gene,a,"\n",sep='\t')                # Show what we have
 
        # Add results to database
      
        ins<-paste0("INSERT INTO ",geneDB," (model,gene,roc,n,p_intercept,p_x,pathway,effect,exp,dataset,mtime) VALUES ('sis~x+stage','",gene,"',",a,",",nrow(d1),",",as.double(as.character(coefficients(s)['(Intercept)','Pr(>|t|)'])),",",          as.double(as.character(coefficients(s)['x','Pr(>|t|)'])),          ",null,",effect,",'",experiment,"','",ld,"',",mt,")")

        con<-dbConnect(MySQL(), user='simon', password='0xford', host='db',dbname='single_cell')
        
        tryCatch(
        {
            executeQuery(ins,db=con)
        },
        error=function(d)
        {
          dbDisconnect(con)
          
          stop(paste("Error",gene,"\n",ins,"\n",sep='\t'))
        },
        finally=dbDisconnect(con))

      },error=function(e)
      {
        cat("Error",gene,"\n",sep='\t')
        print(e);
      },warning=function(e)
      {
        cat("Warning",gene,"\n",sep='\t')
        print(e);
      })
    }
  
    rm(genes)
  }
}

executeQuery("UPDATE ",geneDB," SET pathway=null WHERE pathway='NA'")
