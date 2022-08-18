plot_structure_map_extension <- function(assignments, k, facet, pop_coordinates, sf = NULL, sf_fill_colors = "viridis", sf_line_colors = "viridis",
                                         pop_names = T, viridis.option = "viridis", alt.palette = NULL,
                                         radius_scale = 0.05, label_args = NULL, crop = FALSE,
                                         scale_bar = list(dist = 4, dist_unit = "km", transform = T), compass = list(symbol = 16),
                                         internals){

  
  long <- lat <- pop <- NULL
  
  #===================internals=====================
  for(i in 1:length(internals)){
    assign(names(internals)[i], internals[[i]])
  }
  
  #===================sanity checks=================
  msg <- character()
  pkg.check <- .check.installed("scatterpie")
  pkg.check <- c(pkg.check, .check.installed("sf"))
  if(!is.null(compass) | !is.null(scale_bar)){pkg.check <- c(pkg.check, .check.installed("ggsn"))}
  pkg.check <- c(pkg.check, .check.installed("viridis"))
  
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  
  if(!is.null(sf)){
    use_crs <- sf::st_crs(pop_coordinates)
    
    # polygon palette
    is.poly <- unlist(lapply(sf, function(x) grepl("POLYGON", sf::st_geometry_type(x)[1])))
    poly.sum <- sum(is.poly)
    if(poly.sum > 0){
      if(sf_fill_colors[1] == "viridis"){
        poly_pal <- viridis::viridis(poly.sum, alpha = .2, option = viridis.option)
      }
      else{
        if(length(sf_fill_colors) != poly.sum){
          warning("The length of the provided fill colors is not the same as the number of polygon sf objects to be plotted, defaulting to viridis.\n")
          poly_pal <- viridis::viridis(poly.sum, alpha = .2, option = viridis.option)
        }
        else{
          poly_pal <- sf_fill_colors
        }
      }
      used_poly_pall <- 0
    }
    
    if(sf_line_colors[1] == "viridis"){
      sf_line_colors <- viridis::viridis(length(sf), option = viridis.option)
    }
    else{
      if(length(sf_line_colors) != length(sf)){
        warning("The length of the provided line colors is not the same as the number of sf objects to be plotted, defaulting to viridis.\n")
        sf_line_colors <- viridis::viridis(length(sf), alpha = .2, option = viridis.option)
      }
    }
  }
  
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
  
  # add sf overlay if requested.
  if(!is.null(sf)){
    for(i in 1:length(sf)){
      sf[[i]] <- sf::st_transform(sf[[i]], use_crs)
      
      if(is.poly[i]){
        mp <- mp + ggplot2::geom_sf(data = sf[[i]], fill = poly_pal[used_poly_pall + 1], color = sf_line_colors[i])
        used_poly_pall <- used_poly_pall + 1
      }
      else{
        mp <- mp + ggplot2::geom_sf(data = sf[[i]], color = sf_line_colors[i])
      }
    }
  }
  
  # add the scatterpie
  mp <- mp + ggplot2::theme_bw() +
    scatterpie::geom_scatterpie(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, r = r), cols = colnames(pie_dat)[4:ncol(pie_dat)]) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
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
  
  # scale/compass
  if((!is.null(scale_bar) | !is.null(compass))){
    if(!is.null(sf)){ 
      # make up a null sf object to set extent if needed
      dummy <- sf[[1]]
    }
    else{
      dummy <- pop_coordinates
    }
    
    # grab limit info for cropped data...
    b <- ggplot2::ggplot_build(mp)
    lims <- list(x = b$layout$panel_scales_x[[1]]$limits,
                 y = b$layout$panel_scales_y[[1]]$limits)
    
    
    
    if(!is.null(scale_bar)){
      scale_bar$data <- dummy
      if(crop){
        if(!"anchor" %in% names(scale_bar)){
          scale_bar$anchor <- c(x = lims$x[2], y = lims$y[1] + .1*abs(lims$y[1]))
        }
      }
      mp <- mp + do.call(ggsn::scalebar, args = scale_bar)
    }
    
    
    if(!is.null(compass)){
      compass$data <- dummy
      if(crop){
        if(!"anchor" %in% names(compass)){
          compass$anchor <- c(x = lims$x[2] + abs(lims$x[2]*.05), y = lims$y[2] + abs(lims$y[2]*.05))
        }
      }
      mp <- mp + do.call(ggsn::north, args = compass)
    }
  }
  
  # return
  return(mp)
}
