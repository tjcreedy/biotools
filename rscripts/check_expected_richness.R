#threshold = 0.85
check_expected_richness <- function(map, threshold, method = "chao1", bins = 30){
  require(vegan)
  require(ggplot2)
  require(ggthemes)
  require(cowplot)
  
  stats <- data.frame(cbind(rowSums(map), t(estimateR(map))))
  colnames(stats)[1] <- "reads"
  stats$p_complete <- stats$S.obs/stats[,paste0("S.",method)]
  stats$p_complete_min <- stats$S.obs/(stats[,paste0("S.",method)] + stats[,paste0("se.", method)])
  stats$p_complete_max <- stats$S.obs/(stats[,paste0("S.",method)] - stats[,paste0("se.", method)])
  stats$p_complete_max[stats$p_complete_max > 1] <- 1
  
  passing <- data.frame(x = c(threshold - 0.05, threshold + 0.05),
                        y = c(0.5, 0.5),
                        label = c(paste0(round(sum(stats$p_complete < threshold, na.rm = T)/nrow(stats)*100), "%"),
                                  paste0(round(sum(stats$p_complete > threshold, na.rm = T)/nrow(stats)*100), "%")),
                        col = LETTERS[1:2])
  
  plots1 <- list(ggplot(stats, aes(x = p_complete))+
                  labs(y = "Count of samples") +
                  geom_histogram(bins = bins),
                ggplot(stats)+
                  labs(y = "Cumulative proportion of samples <= x")+
                  stat_ecdf(aes(x = p_complete))+
                  stat_ecdf(aes(x = p_complete_min), linetype = 2) + 
                  stat_ecdf(aes(x = p_complete_max), linetype = 2) + 
                  geom_text(data = passing, aes(x = x, y = y, label = label, col = col))+
                  scale_colour_manual(values = c("red", "green"), guide = F))
  plots2 <- list(ggplot(stats, aes(x = reads, y = p_complete, ymin = p_complete_min, ymax = p_complete_max))+
                   labs(x = "Number of reads") +
                   scale_x_continuous(trans = 'sqrt') +
                   geom_errorbar(alpha = .3) +
                   geom_point(alpha = .6),
                 ggplot(stats, aes(x = S.obs, y = p_complete, ymin = p_complete_min, ymax = p_complete_max))+
                   labs(x = "Observed richness") +
                   scale_x_continuous(trans = 'sqrt') +
                   geom_errorbar(alpha = .3) +
                   geom_point(alpha = .6))
  
  plots <- c(lapply(plots1, function(p){
    p+
      labs(x = "Proportion of expected OTUs observed")+
      geom_vline(xintercept = threshold, col = "red")+
      theme_tufte()
  }),
    lapply(plots2, function(p){
    p+
      labs(y = "Proportion of expected OTUs observed")+
      geom_hline(yintercept = threshold, col = "red")+
      theme_tufte()
  }))
  
  print(plot_grid(plots[[1]]+
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank()), 
                  plots[[3]]+
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank()),
                  plots[[2]],
                  plots[[4]],
                  nrow = 2, align = "v")
  )
  
  output = list(
    minpass = row.names(stats[stats$p_complete_min >= threshold,]),
    meanpass = row.names(stats[stats$p_complete >= threshold,]),
    maxpass = row.names(stats[stats$p_complete_max >= threshold,])
  )
  
  return(list(passnames = output, stats = stats, plots = plots))
}