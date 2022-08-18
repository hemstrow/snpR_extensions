
# snpR internals packaged here for use:
.check.installed <- function(pkg, install.type = "basic", source = NULL){
  if(!pkg %in% rownames(utils::installed.packages())){
    say <- paste0("Package '", pkg, "' not found.")
    cat(say, "")
    cat("Install? (y or n)\n")
    
    resp <- readLines(n = 1)
    resp <- tolower(resp)
    
    
    if(resp != "y"){
      stop(say)
    }
    else{
      if(install.type == "basic"){
        utils::install.packages(pkg)
      }
      if(install.type == "bioconductor"){
        if(!"BiocManager" %in% utils::installed.packages()){
          utils::install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      }
      if(install.type == "github"){
        if(!"remotes" %in% utils::installed.packages()){
          utils::install.packages("remotes")
        }
        remotes::install_github(source)
      }
      return(TRUE)
    }
  }
  else{return(TRUE)}
}


.check.snpR.facet.request <- function(x, facets, remove.type = "snp", return.type = FALSE, fill_with_base = TRUE, return_base_when_empty = TRUE){
  if(any(facets == "all")){
    facets <- x@facets
  }
  if(is.null(facets) | isFALSE(facets)){
    facets <- ".base"
    if(return.type){
      return(list(facets, ".base"))
    }
    else{
      return(".base")
    }
  }
  
  
  # remove the facet parts as requested.
  facets <- .split.facet(facets)
  to.remove <- logical(length(facets))
  missing.facets <- character(0)
  facet.types <- character(length(facets))
  for(i in 1:length(facets)){
    if(identical(facets[[i]], ".base")){next()}
    facets[[i]] <- sort(facets[[i]])
    samp.facets <- which(facets[[i]] %in% colnames(x@sample.meta))
    snp.facets <- which(facets[[i]] %in% colnames(x@snp.meta))
    missing.facets <- c(missing.facets, facets[[i]][which(!((1:length(facets[[i]])) %in% c(samp.facets, snp.facets)))])
    
    if(remove.type == "snp"){
      if(length(snp.facets) > 0){
        if(length(snp.facets) == length(facets[[i]])){
          to.remove[i] <- T
        }
        else{
          facets[[i]] <- facets[[i]][-snp.facets]
        }
      }
    }
    
    else if(remove.type == "sample"){
      if(length(samp.facets) > 0){
        if(length(samp.facets) == length(facets[[i]])){
          to.remove[i] <- T
        }
        else{
          facets[[i]] <- facets[[i]][-samp.facets]
        }
      }
    }
    
    else if(remove.type == "complex"){
      if(length(snp.facets) > 0 & length(samp.facets) > 0){
        to.remove[i] <- T
      }
    }
    
    else if(remove.type == "simple"){
      if(!(length(snp.facets) > 0 & length(samp.facets) > 0)){
        to.remove[i] <- T
      }
    }
    
    if(return.type){
      samp.facets <- which(facets[[i]] %in% colnames(x@sample.meta))
      snp.facets <- which(facets[[i]] %in% colnames(x@snp.meta))
      missing.facets <- c(missing.facets, facets[[i]][which(!((1:length(facets[[i]])) %in% c(samp.facets, snp.facets)))])
      
      if(length(samp.facets) > 0 & length(snp.facets) > 0){
        facet.types[i] <- "complex"
      }
      else if(length(samp.facets) > 0){
        facet.types[i] <- "sample"
      }
      else if(length(snp.facets) > 0){
        facet.types[i] <- "snp"
      }
    }
    facets[[i]] <- paste(facets[[i]], collapse = ".")
  }
  
  if(length(missing.facets) > 0){
    dups <- which(duplicated(missing.facets))
    if(length(dups) > 0){
      stop("Facet(s) ", paste(missing.facets[-dups], collapse = ", "), " not found in x metadata.\n")
    }
    stop("Facet(s) ", paste(missing.facets, collapse = ", "), " not found in x metadata.\n")
  }
  
  if(any(to.remove)){
    facets <- facets[-which(to.remove)]
    facet.types <- facet.types[-which(to.remove)]
    if(fill_with_base){
      if(!".base" %in% facets){
        facets <- c(facets, ".base")
        facet.types <- c(facet.types, ".base")
      }
    }
  }
  
  # fix if we've removed everything, return the base facet if return_base_when_empty is TRUE
  if(return_base_when_empty){
    if(length(facets) == 0){
      facets <- ".base"
      if(return.type){
        return(list(facets, ".base"))
      }
      else{
        return(".base")
      }
    }
  }
  
  
  # fix the facet type for .base
  if(return.type){
    base.facets <- which(facet.types == "")
    if(length(base.facets) > 0){
      facet.types[base.facets] <- ".base"
    }
  }
  
  # remove duplicates and return
  dups <- duplicated(facets)
  if(any(dups)){
    facets <- facets[-which(dups)]
    facet.types <- facet.types[-which(dups)]
  }
  if(return.type){
    return(list(unlist(facets), facet.types))
  }
  else{
    return(unlist(facets))
    
  }
}

.split.facet <- function(facet) strsplit(facet, "(?<!^)\\.", perl = T)

.add.facets.snpR.data <- function(x, facets = NULL){
  #===========================binding issues for unquoted variables============
  .snp.id <- facet <- subfacet <- NULL
  
  
  if(is.null(facets[1])){return(x)}
  facets <- .check.snpR.facet.request(x, facets)
  if(is.null(facets[1])){return(x)}
  #===========================turn into list========
  # need to fix any multivariate facets (those with a .)
  comp.facets <- grep("(?<!^)\\.", facets, perl = T)
  if(length(comp.facets) != 0){
    run.facets <- as.list(facets[-c(comp.facets)])
    facet.list <- c(run.facets, .split.facet(facets[comp.facets]))
  }
  else{
    facet.list <- as.list(facets)
  }
  
  #===========================sanity checks=========
  # check that we haven't done these facets before, and remove any that we have.
  all.facets <- character(length = length(facet.list))
  for(i in 1:length(facet.list)){
    all.facets[i] <- paste0(facet.list[[i]], collapse = ".")
  }
  if(all(all.facets %in% x@facets)){return(x)}
  else if (any(all.facets %in% x@facets)){
    already.added <- all.facets[which(all.facets %in% x@facets)]
    facet.list <- facet.list[-which(all.facets %in% x@facets)]
  }
  
  # check that all of the facet columns actually exist
  all.facets <- unlist(unlist(facet.list))
  opts <- c(colnames(x@sample.meta), colnames(x@snp.meta))
  used <- opts[which(opts %in% all.facets)]
  if(!all(all.facets %in% opts)){
    bad.facets <- which(!(all.facets %in% c(colnames(x@sample.meta), colnames(x@snp.meta))))
    stop(paste0("Facets ", paste0(all.facets[bad.facets], collapse = " "), " not found in sample or snp metadata.\n"))
  }
  if(any(duplicated(used))){
    stop(paste0("Facets ", paste0(used[which(duplicated(used))], collapse = " "), " are duplicated in the sample and or snp metadata.\n"))
  }
  
  #===========================process each facet.===================
  # for sample acets, prep summary data and geno tables.
  # For snp facets, nothing to do here besides add them to the facet list and facet type list.
  # For complex facets, prep summary tables for the sample facet portion if they don't exist and add to facet list.
  added.facets <- character(0)
  oac <- cbind(data.table::as.data.table(x@facet.meta[,c("facet", "subfacet", ".snp.id")]),
               data.table::as.data.table(x@ac)) # grab original ac with meta for later
  for(k in 1:length(facet.list)){
    facets <- facet.list[[k]] # column levels for this facet.
    #=========================figure out unique levels for the facet==========
    # figure out what kind of facets we are working with.
    
    # save info and get the unique options for each facet
    # save info
    if(length(facets) > 1){
      x@facets <- c(x@facets, paste0(facets, collapse = "."))
      added.facets <- c(added.facets, paste0(facets, collapse = "."))
    }
    else{
      x@facets <- c(x@facets, facets)
      added.facets <- c(added.facets, facets)
    }
    x@facet.type <- c(x@facet.type, "sample")
    
    
    # get unique options for this facet
    sample.meta <- x@sample.meta[colnames(x@sample.meta) %in% facets]
    sample.meta <- sample.meta[,sort(colnames(sample.meta))]
    if(!is.data.frame(sample.meta)){
      sample.meta <- as.data.frame(sample.meta, stringsAsFactors = F)
      colnames(sample.meta) <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% facets]
    }
    sample.meta <- dplyr::mutate_all(sample.meta, as.character) # this fixes some really obscure bugs with integers in columns.
    sample.opts <- unique(sample.meta)
    if(!is.data.frame(sample.opts)){
      sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
      colnames(sample.opts) <- facets[which(facets %in% colnames(x@sample.meta))]
    }
    sample.opts <- dplyr::arrange_all(sample.opts)
    
    
    gs <- x@geno.tables
    #=========================get gs matrices==========
    for(i in 1:nrow(sample.opts)){
      matches <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[i,]))))
      t.x <- genotypes(x)[,matches]
      tgs <- .tabulate_genotypes(t.x, x@mDat)
      gs$gs <- plyr::rbind.fill.matrix(gs$gs, tgs$gs)
      gs$as <- plyr::rbind.fill.matrix(gs$as, tgs$as)
      gs$wm <- plyr::rbind.fill.matrix(gs$wm, tgs$wm)
      # fix NAs that show up when there are less called genotype options in one facet level than in all levels!
      if(ncol(gs$gs) != ncol(tgs$gs)){
        gs <- lapply(gs, function(x){x[is.na(x)] <- 0;x})
      }
      
      x@facet.meta <- rbind(x@facet.meta,
                            cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(tgs$gs)),
                                             subfacet = rep(paste0(sample.opts[i,], collapse = "."), nrow(tgs$gs)),
                                             facet.type = rep("sample", nrow(tgs$gs)), stringsAsFactors = F),
                                  x@snp.meta))
    }
    
    
    #=========================sort, pack, and return==========
    # sort
    x@facet.meta <- dplyr::mutate_if(.tbl = x@facet.meta, is.factor, as.character)
    x@facet.meta$.reorder <- 1:nrow(x@facet.meta)
    x@facet.meta <- dplyr::arrange(x@facet.meta, .snp.id, facet, subfacet)
    gs$gs <- gs$gs[x@facet.meta$.reorder,, drop = F]
    gs$as <- gs$as[x@facet.meta$.reorder,, drop = F]
    gs$wm <- gs$wm[x@facet.meta$.reorder,, drop = F]
    x@facet.meta <- x@facet.meta[,-ncol(x@facet.meta)]
    
    # output
    x@geno.tables <- gs
  }
  # add and sort ac formated data.
  .make_it_quiet(nac <- format_snps(x, output = "ac", facets = added.facets))
  nac <- data.table::as.data.table(nac)
  nac <- rbind(oac, nac[,c("facet", "subfacet", ".snp.id", "n_total","n_alleles", "ni1", "ni2")])
  nac <- dplyr::mutate_if(.tbl = nac, is.factor, as.character)
  nac <- dplyr::arrange(nac, .snp.id, facet, subfacet)
  nac <- as.data.frame(nac)
  x@ac <- nac[,-c(1:3)]
  
  # add in dummy rows to stats
  sm <- data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% added.facets, c("facet", "subfacet", "facet.type", colnames(x@snp.meta))])
  if(ncol(x@stats) - ncol(sm) > 0){
    sm <- cbind(sm, matrix(NA, nrow(sm), ncol(x@stats) - ncol(sm)))
  }
  colnames(sm) <- colnames(x@stats)
  os <- data.table::as.data.table(x@stats)
  if(ncol(os) - ncol(sm) > 0){
    os <- rbind(os, cbind(sm, matrix(NA, nrow(sm), ncol(os) - ncol(sm))))
  }
  else{
    os <- rbind(os, sm)
  }
  
  os <- dplyr::mutate_if(.tbl = os, is.factor, as.character)
  x@stats <- as.data.table(dplyr::arrange(os, .snp.id, facet, subfacet))
  
  
  return(x)
}



.tabulate_genotypes <- function(x, mDat, verbose = F){
  
  # fix for if x is a vector (only one individual) and convert to data.table
  if(!is.data.frame(x)){
    x <- data.frame(samp = x)
  }
  x <- data.table::setDT(x)
  
  
  # get a genotype table
  snp_form <- nchar(x[1,1])   # get information on data format
  x <- data.table::melt(data.table::transpose(x, keep.names = "samp"), id.vars = "samp") # transpose and melt
  
  gmat <- data.table::dcast(data.table::setDT(x), variable ~ value, value.var='value', length) # cast
  gmat <- gmat[,-1]
  mis.cols <- -which(colnames(gmat) == mDat)
  if(length(mis.cols) > 0){
    tmat <- gmat[,mis.cols, with = FALSE] # remove missing data
  }
  else{
    tmat <- gmat
  }
  
  #get matrix of allele counts
  #initialize
  hs <- substr(colnames(tmat),1,snp_form/2) != substr(colnames(tmat), (snp_form/2 + 1), snp_form*2) # identify heterozygotes.
  if(verbose){cat("Getting allele table...\n")}
  as <- unique(unlist(strsplit(paste0(colnames(tmat)), "")))
  amat <- data.table::as.data.table(matrix(0, nrow(gmat), length(as)))
  colnames(amat) <- as
  
  #fill in
  for(i in 1:length(as)){
    b <- grep(as[i], colnames(tmat))
    hom <- which(colnames(tmat) == paste0(as[i], as[i]))
    if(length(hom) == 0){
      het <- b
      set(amat, j = i, value = rowSums(tmat[,het, with = FALSE]))
    }
    else{
      het <- b[b != hom]
      if(length(het) > 0){
        if(data.table::is.data.table(tmat[,het, with = FALSE])){
          set(amat, j = i, value = (tmat[,hom, with = FALSE] * 2) + rowSums(tmat[,het, with = FALSE]))
        }
        else{
          amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
        }
      }
      else{
        set(amat, j = i, value = (tmat[,hom, with = FALSE] * 2))
      }
    }
  }
  return(list(gs = as.matrix(tmat), as = as.matrix(amat), wm = as.matrix(gmat)))
}


.make_it_quiet <- function(fun){
  return(invisible(utils::capture.output(fun)))
}