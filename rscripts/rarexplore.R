rarexplore <- function(map, size, seed = 12345){
  require(ggplot2)
  require(reshape2)
  require(plyr)
  require(magrittr)
  require(scales)
  require(vegan)
  raresample = size[2]
  rared <- lapply(size, function(raresample){
    set.seed(seed)
    raresample <- round(raresample)
    rared <- suppressWarnings(rrarefy(map, raresample))
    passrows <- rowSums(rared) >= raresample
    passcols <- colSums(rared) > 0
    rared <- rared[passrows, passcols] %>%
      matrix(., nrow = sum(passrows), ncol = sum(passcols),
             dimnames = list(names(passrows)[passrows], names(passcols)[passcols]))
    return(setNames(c(raresample,dim(rared)),
                    c("raresamples","Samples","OTUs"))) 
  }) %>% do.call(rbind, .) %>%
    data.frame() %>%
    melt(., id.vars = "raresamples", variable.name = "measure", value.name = "value")
  
  rared$value <- ifelse(rared$measure == "Samples", rared$value/nrow(map), rared$value/ncol(map))
  
  ggplot(data = rared, aes(x = raresamples, y = value))+
    geom_point()+
    theme_bw()+ 
    geom_line()+ 
    scale_y_continuous(labels = scales::percent)+
    labs(x = "Rarefaction sample", y = "Percentage of total remaining")+
    facet_wrap(~measure, scales = "free_y")
    
}
