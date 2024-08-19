options(error=traceback) # Switch on debug if error

# Parse command line args

args = commandArgs(trailingOnly=TRUE) # Accept methods to try through Rscript command line

batchToDo<-as.integer(args[1])
if(!is.na(batchToDo)) args<-args[-1]

verbosity<-F

if(!is.na(args[1])&&(args[1]=='-v'))
{
  verbosity<-T
  args<-args[-1]
  
  message("*** VERBOSE MODE ***")
}
  
# Load libraries, get data & set seed for reproducibility ---------------------

set.seed(123)    # seed for reproducibility

require(dplyr,quietly = TRUE,warn.conflicts=FALSE)   # for data cleaning
require(pROC,quietly = TRUE,warn.conflicts=FALSE)
require(caret,quietly = TRUE,warn.conflicts=FALSE)
require(VGAM,quietly = TRUE,warn.conflicts=FALSE)
require(R.utils,quietly = TRUE,warn.conflicts=FALSE)

if(file.exists("sql.r")) source("sql.r") else source('~/R/sql.r')

source('znrf3.r')

source('ml_tools.r')

# Set up parallel environment

require(doParallel,quietly=TRUE,warn.conflicts=FALSE)

#

mentalMethods<-c('neuralnet','bagEarthGCV','gam','xgbLinear','xgbTree','svmRadialCost','extraTrees','xgbDART','widekernelpls','rqnc','nodeHarvest','SBC','gaussprRadial','DENFIS','qrnn','gamLoess','knn')

# Set training control
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = verbosity)

`%ni%` <- Negate(`%in%`)

# Low hanging fruit-a-thon

todos<-c(60)
if(!is.na(batchToDo)) todos<-c(batchToDo)

for(to in todos) # Go for low hanging fruit first
{
  for(g in c('Sertoli_X_Y',"Granulosa_X_Y")) # Model for each cell type
  {
    # Load data and respective SIS scores
    
    if(file.exists(paste0('sex_gene_all_Normal_all_',g,'.rds'))) cellset<-readRDS(paste0('sex_gene_all_Normal_all_',g,'.rds')) else cellset<-readRDS(paste0(modelDir,'/sex_gene_all_Normal_all_',g,'.rds'));
  
    # Trim the cellset to the appropriate wild-types

    methods<-getMethodsToDo(g,to)
    if(length(args)>0) methods<-args # methods[args %in% getMethodsToDo(g,to)]
    
    if(g=='Sertoli_X_Y') 
    {
      cellset<-cellset[cellset$genotype=='xy',] # Only xy wt data
      sis<-ifelse(cellset[,'cell_type_s'],1,0)  # Use single genotype cell type
      training<-cellset$training_s
      todo<-methods
    } else
    {
      cellset<-cellset[cellset$genotype=='xx',] # Only xx wt data
      sis<-ifelse(cellset[,'cell_type_g'],1,0)  # Use single genetype cell type
      training<-cellset$training_g
      
      todo<-methods
    }
    
    # Get data ready for caret
    
    to.remove <- c("x","ensemble","cell","genotype","label","indiv_labels","cell_type","cell_type_s","cell_type_g","training_s","training_g","gene","dev_gene","house_keeping","d$cell")
    x<-subset(cellset,select = names(cellset) %ni% to.remove)
    
    ##################################
    # Loop through available methods #
    ##################################

    # Either do genotype models with no timeout >= to and no 'ok' record at any time (it's done)
    # If we supply methods via command line, do only those - if needed
  
    if(file.exists('/NGS/working_projects/scRNA-Seq/ML_Results_20201212/')) methodDir<-'/NGS/working_projects/scRNA-Seq/ML_Results_20201212/methods/' 

    for(method in methods)
    {
      if(!method %in% todo) next

      if(method %in% veto) 
      {
         message(paste0("Method ",method," is in the veto list"))
	 next
      }


      cores<-ifelse(method %in% mentalMethods,1,parallel::detectCores())
      
      registerDoParallel(cores=cores)
      
      message(paste("Using ",cores," cores."))
      
      modelFile <-paste0(methodDir,g,'_',method,'.rds'); # Model filed
      
      message(paste("Looking at ",method,g,to,sep='\t'))

      dbok<-singleQuery("SELECT COUNT(*) FROM ",methodDB," WHERE method='",method,"' AND genotype='",g,"' AND error='ok' AND roc_auc IS NULL")>0
      modelok<-file.exists(modelFile)
      
      tryCatch( # Might get an error, so catch it to avoid total crash
      {
	      message(paste0(modelFile,"\t",file.exists(modelFile)))

        if(!modelok)
        {
          # Train the model
        
          message(paste("\tRunning ",method,g,to,sep='\t'))
          
          #withTimeout({    
            if(verbosity)
            {
              model <- train(sis ~ .,
                        data = cbind(sis, x),
                        method = method,
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        subset=training,
                        trControl = train_control)
              
            } else
            {
              capture.output(suppressMessages(suppressWarnings(model <- train(sis ~ .,
                          data = cbind(sis, x),
                          method = method,
                          preProcess = c("center", "scale"),
                          tuneLength = 10,
                          subset=training,
                          trControl = train_control))),file='/dev/null')
            }
          #},timeout=to,cpu=9999999,onTimeout="error")  
          
          message(paste0("\tModel run completed - saving model to '",modelFile,"'"))
            
          ssaveRDS(model,file=modelFile) # Save the model
          
        } else
	      {
		      message("\tLoading model");

          model<-readRDS(modelFile);
	      }

        # Update database
        
        if(!is.na(model)&&!dbok)
        {
	        message("\tUpdating database")

          sis_p <- predict(model, x) # Get predicted values
          rsq <- cor(sis,sis_p)^2    # R2
          
          r<-roc(cellset$sis>0.5,sis_p,quiet=TRUE); # TODO: Split data?
          a<-auc(r)    # ROC AUC
          
          if(to>0)
          {
            db<-connect()
            sql<-paste0("REPLACE INTO ",methodDB," (method,timeout,genotype,error,r2,roc_auc) VALUES ('",method,"',",to,",'",g,"','ok',",ifelse(is.na(rsq),"null",rsq),",",ifelse(is.na(a),"null",a),")")
            #message(sql)
            dbExecute(conn=db,statement=sql)
            dbDisconnect(db)            
          }
          
          message(paste("\t",g,method,rsq,a,sep='\t'))
        }
      }
      ,error=function(e)
      {
        saveRDS(e,file ='/scratch/e.rds')
        
        message(paste("\tError",method,e))
        
        if(to>0)
        {
          db<-connect()
          sql<-paste0("REPLACE INTO ",methodDB," (method,timeout,genotype,error) VALUES ('",method,"',",to,",'",g,"','",escape(as.character(e)),"')")
          #message(sql)
          dbExecute(conn=db,statement=sql)
          dbDisconnect(db)            
        }
      })
    }
  }
}
