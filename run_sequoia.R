run_sequoia_extension <- function(x, facets = NULL, run_dupcheck = FALSE, 
                                  run_parents = FALSE, run_pedigree = FALSE, 
                                  min_maf = 0.3, min_ind = 0.5, ...,
                                  internals){
  
  #===================internals=====================
  for(i in 1:length(internals)){
    assign(names(internals)[i], internals[[i]])
  }
  
  #===================sanity checks===============
  # check that provided snpRdata objects are in the correct format
  if(!snpR::is.snpRdata(x)){
    stop("Not a snpRdata object.\n")
  }
  
  
  msg <- character(0)
  warn <- character(0)
  
  if(run_pedigree & !run_parents){
    msg <- c(msg, "Parents must be run before pedigree construction!\n")
  }
  
  if(min_maf < 0.25){
    warn <- c(warn, "Minor allele frequencies below 0.25 are not recommended for sequoia.\n")
  }
  if(min_ind < 0.5){
    warn <- c(warn, "Genotypes sequenced in fewer than 50% of individuals will automatically be removed by Sequoia.\n")
  }
  
  req.columns <- c("ID", "Sex")
  if(any(!req.columns %in% colnames(sample.meta(x)))){
    msg <- c(msg, "Columns named 'ID' and 'Sex' are required in the sample metadata.\n")
  }
  
  # check columns in metadata
  col_logi <- sum(c("BirthYear" %in% colnames(sample.meta(x)),
                    sum(c("BYmin", "BYmax") %in% colnames(sample.meta(x))) == 2))
  
  if(col_logi == 0){
    msg <- c(msg, "Columns named either BirthYear or both BYmin and BYmax are required in sample meta.\n")
  }
  else if(col_logi == 2){
    warn <- c("Both BirthYear, BYmin, and BYmax are contained in the sample metadata, defaulting to BirthYear. To change this behavior, remove the BirthYear column.\n")
  }
  
  if(length(warn) > 0){
    warning(warn)
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #===================prep=========================
  # prep facets
  facets <- .check.snpR.facet.request(x, facets, "none")
  x <- .add.facets.snpR.data(x, facets)
  tasks <- .get.task.list(x, facets = facets)
  
  # initialize
  out <- vector("list", nrow(tasks))
  
  #==================run===========================
  for(i in 1:nrow(tasks)){
    
    names(out)[i] <- paste0(tasks[i,], collapse = "_")
    cat("Running facet: ", paste0(tasks[i,c(1,3)], collapse = " "), "\tsubfacet: ", paste0(tasks[i,c(2,4)], collapse = " "))
    
    
    # subset the data
    tmatches <- .fetch.sample.meta.matching.task.list(x, tasks[i,])
    if(tasks[i,3] == ".base"){
      snpmatches <- 1:nsnps(x)
    }
    else{
      snpmatches <- which(snp.meta(x)[,tasks[i, 3]] == tasks[i,4])
    }
    tdat <- x[snpmatches, tmatches]
    
    
    
    # filter snps for running sequoia - requires high MAF >0.25 and 50% inds otherwise inds dropped from seq
    tdat <- filter_snps(x = tdat, maf = min_maf, min_ind = min_ind, re_run = "full") # set as arguments later
    
    # format the snps for sequoia
    tdat <- format_snps(x = tdat, output = "sequoia")
    
    # initialize output
    out[[i]] <- vector("list")
    
    # run sequoia
    if(run_dupcheck){
      dups <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, Module = "dup")
      out[[i]]$dups <- dups
    }
    
    if(run_parents){
      out[[i]]$parents <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, Module = "par", ...)
    }
    
    if(run_pedigree){
      out[[i]]$pedigree <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, SeqList = out[[i]]$parents, Module = "ped", ...) 
    }
  }
  return(out)
}
