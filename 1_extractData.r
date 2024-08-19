require(dplyr,quietly = TRUE,warn.conflicts=FALSE)

source('znrf3.r')

# Takes 30mins

memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))

if(memfree<32000000)
{
  stop("Not enough memory - use a bigger server!")
}

source('normalize_expr_data.R');

# Load Richard's data

is_sparse_matrix <- function(x)
{
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}

if(!exists('cdsg')) cdsg<-readRDS(rawData2[1]);
if(!exists('cdss')) cdss<-readRDS(rawData2[2]);

for(rd in c(rawData,rawData2))
{
  ld<-add(rd,"_long")
  
  if(!file.exists(rd))
  {
    message("ERROR: File not found ",rd)
    next
  }
  
  message(paste(rd,ld))
  
  #if(older(rd,ld))
  {
    message(paste0("Loading data from ",rd));
    
    cds <- readRDS(rd)

    cts<-c()

    g<-cellify(colnames(cdsg[,cdsg$cell_type=='Granulosa']),cdsg[,cdsg$cell_type=='Granulosa']$s_name)
    s<-cellify(colnames(cdss[,cdss$cell_type=='Sertoli']),cdss[,cdss$cell_type=='Sertoli']$s_name)
    c<-cellify(colnames(cds),cds$s_name)
     
    cds$cell_type_g<-c %in% g
    cds$cell_type_s<-c %in% s

    gnotfound<-sum(cellify(colnames(cdsg[,cdsg$cell_type=='Granulosa']),cdsg[,cdsg$cell_type=='Granulosa']$s_name) %notin% cellify(colnames(cds),cds$s_name))
    snotfound<-sum(cellify(colnames(cdss[,cdsg$cell_type=='Sertoli']),cdss[,cdss$cell_type=='Sertoli']$s_name) %notin% cellify(colnames(cds),cds$s_name))
    
    if(rd==rawData) message(paste("\tCells not found in combi:",gnotfound,snotfound,sep='\t'))
    
    cds$training_g<-cellify(colnames(cds),cds$s_name) %in% cellify(colnames(cdsg),cdsg$s_name)
    cds$training_s<-cellify(colnames(cds),cds$s_name) %in% cellify(colnames(cdss),cdss$s_name)
    
    message("\tNormalising");

    nm <- normalize_expr_data(cds)

    message("\tCreating long format data frame")
  
    d<-as.data.frame(summary(nm));
    
    nm2<-summary(normalize_expr_data(cds,norm_method='none'))
    nm2$raw<-nm2$x
    
    d$raw<-nm2[,'raw']

    d$ensemble<-nm@Dimnames[[1]][d$i]
    d$cell<-nm@Dimnames[[2]][d$j]
    #d$sample_index<-as.integer(substring(d$cell,20))
    m<-match(d$cell,rownames(cds@colData))
  
    d$genotype<-cds@colData[m,'genotype']
    d$stage<-cds@colData[m,'stage']  
    d$label<-cds@colData[m,'label']
    d$indiv_labels<-cds@colData[m,'indiv_label']  
    #d$somite<-cds@colData[m,'somite']
    d$cell_type<-cds@colData[m,'cell_type']
    d$cell_type_s<-cds@colData[m,'cell_type_s']
    d$cell_type_g<-cds@colData[m,'cell_type_g'] 
    d$training_s<-cds@colData[m,'training_s']
    d$training_g<-cds@colData[m,'training_g'] 
  
    genes<-fData(cds)
    m<-match(d$ensembl,rownames(genes))
    d$gene<-genes[m,'gene_short_name']
  
    rm(nm,cds,m,genes); gc(verbose=FALSE);
  
    d<-d[ , -which(names(d) %in% c("i","j"))]
  
    d$ensemble<-as.factor(d$ensemble)
    d$cell<-as.factor(d$cell)
    d$genotype<-as.factor(d$genotype)
    d$stage<-as.factor(d$stage)
    d$label<-as.factor(d$label)
    d$indiv_labels<-as.factor(d$indiv_labels)
    #d$somite<-as.factor(d$somite)
    d$cell_type<-as.factor(d$cell_type)
    d$gene<-as.factor(d$gene)
  
    # Mark known develomental genes
  
    d$dev_gene<-'no';
    d[d$gene %in% ovaryDevGenes,'dev_gene']='ovary'
    d[d$gene %in% testesDevGenes,'dev_gene']='testes'
  
      # Mark housekeeping genes
  
    load(houseKeepingData);
    d$house_keeping<-d$gene %in% Mouse_HK_genes$'Gene'
  
    gc();
 
    stop("ready to save") 
    message("\tSaving long data")

    saveRDS(file=ld,object=d)
  } else
  {
    message(paste0("Already present and newer:\t",ld))
  }
}
