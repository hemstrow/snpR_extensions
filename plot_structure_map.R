plot_structure_map_extension <- function(assignments, k, facet, pop_coordinates, layers = NULL,
                                         pop_names = T, viridis.option = "viridis", alt.palette = NULL,
                                         radius_scale = 0.05, label_args = NULL, crop = FALSE,
                                         scale_bar = list(), 
                                         compass = list(style = ggspatial::north_arrow_fancy_orienteering(), location = "br"),
                                         internals){

  
  long <- lat <- pop <- NULL
  
  #===================internals=====================
  # for(i in 1:length(internals)){
  #   assign(names(internals)[i], internals[[i]])
  # }
  
  #===================sanity checks=================
  msg <- character()
  use_crs <- sf::st_crs(pop_coordinates)
  
  K_opts <- unique(assignments$plot_data$K)
  K_opts <- as.numeric(gsub("K = ", "", K_opts))
  if(!k %in% K_opts){
    msg <- c(msg, "Requested value of k not found in provided assignments.\n")
  }
  
  if(length(msg) > 0){stop(msg, "\n")}
  
  #==================prep for plot ===========================
  # generate plotting data.frame
  pie_dat <- as.data.frame(matrix(0, nrow = length(unique(assignments$plot_data[,which(colnames(assignments$plot_data) == facet)])), ncol = 3 + k))
  colnames(pie_dat) <- c("pop", "lat", "long", paste0("Cluster ", 1:k))
  tpd <- assignments$plot_data[assignments$plot_data$K == paste0("K = ", k),]
  tpd <- tpd[,c(facet, "Cluster", "Percentage")]
  tpd$Cluster <- as.numeric(tpd$Cluster)
  anc <- tapply(tpd$Percentage, tpd[,c(facet, "Cluster")], mean)
  
  ## get the pie coordinates
  if(nrow(pop_coordinates) != nrow(pie_dat)){
    stop(paste0("The number of unique options of ", facet, " (", nrow(pie_dat), ") is not equal to the number of provided coordinates (", nrow(pop_coordinates), ").\n"))
  }
  pie_dat[,1] <- as.data.frame(pop_coordinates)[,facet]
  pie_dat[,3:2] <- sf::st_coordinates(pop_coordinates)
  pie_dat[,4:ncol(pie_dat)] <- anc[match(pie_dat[,1], row.names(anc)),]
  
  # figure out the radius to use
  lat_range <- range(pie_dat$lat)
  long_range <- range(pie_dat$long)
  r <- min(lat_range[2] - lat_range[1], long_range[2] - long_range[1])
  r <- r*radius_scale
  
  
  #============make the plot====================
  mp <- ggplot2::ggplot()
  
  # pre-add other layers if requested
  if(!is.null(layers)){
    for(i in 1:length(layers)){
      mp <- mp + layers[[i]]
    }
  }
  
  # add the scatterpie
  mp <- mp + ggplot2::theme_bw() +
    scatterpie::geom_scatterpie(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, r = r), cols = colnames(pie_dat)[4:ncol(pie_dat)]) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Ancestry Cluster")) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  if(crop){
    xr <- c(min(pie_dat$long - r), max(pie_dat$long) + r)
    yr <- c(min(pie_dat$lat - r), max(pie_dat$lat) + r)
    mp <- mp +
      ggplot2::xlim(xr) +
      ggplot2::ylim(yr)
  }
  
  if(!is.null(alt.palette)){
    mp <- mp + ggplot2::scale_fill_manual(values = alt.palette)
  }
  else{
    mp <- mp + ggplot2::scale_fill_viridis_d(option = viridis.option)
  }
  if(pop_names){
    # add labels
    #mp + ggplot2::geom_label(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, label = pop), size = r*label_scale)
    label_call <- c(list(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, label = pop)), label_args)
    mp <- mp + do.call(ggrepel::geom_label_repel, label_call)
  }
  
  if(is.list(scale_bar)){
    if(length(scale_bar) > 0){
      mp <- mp + do.call(ggspatial::annotation_scale, args = scale_bar)
    }
    else{
      mp <- mp + ggspatial::annotation_scale()
    }
    
  }
  if(is.list(compass)){
    if(length(compass) > 0){
      mp <- mp + do.call(ggspatial::annotation_north_arrow, args = compass)
    }
    else{
      mp <- mp + ggspatial::annotation_north_arrow()
    }
  }
  
  # return
  return(mp)
}
