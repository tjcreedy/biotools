#threshold = 0.85
check_expected_richness <- function(map, threshold, method = "chao1", bins = 30){
  require(vegan)
  require(ggplot2)
  require(ggthemes)
  require(cowplot)
  
  stats <- data.frame(cbind(rowSums(map), t(estimateR(map))))
  colnames(stats)[1] <- "reads"
  stats$p_complete <- stats$S.obs/stats[,paste0("S.",method)]
  
  passing <- data.frame(x = c(threshold - 0.05, threshold + 0.05),
                        y = c(0.5, 0.5),
                        label = c(paste0(round(sum(stats$p_complete < threshold, na.rm = T)/nrow(stats)*100), "%"),
                                  paste0(round(sum(stats$p_complete > threshold, na.rm = T)/nrow(stats)*100), "%")),
                        col = LETTERS[1:2])
  
  plots <- list(ggplot(stats, aes(x = p_complete))+
                  labs(y = "Count of samples")+
                  geom_histogram(bins = bins),
                ggplot(stats, aes(x = p_complete))+
                  labs(y = "Cumulative proportion of samples <= x")+
                  stat_ecdf()+
                  geom_text(data = passing, aes(x = x, y = y, label = label, col = col))+
                  scale_colour_manual(values = c("red", "green"), guide = F),
                ggplot(stats, aes(x = p_complete, y = reads))+
                  labs(y = "Number of reads")+
                  geom_point(),
                ggplot(stats, aes(x = p_complete, y = S.obs))+
                  labs(y = "Number of OTUs")+
                  geom_point())
  
  plots <- lapply(plots, function(p){
    p+
      labs(x = "Proportion of expected OTUs observed")+
      geom_vline(xintercept = threshold, col = "red")+
      theme_tufte()
  })
  
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
  
  return(list(row.names(stats[stats$p_complete >= threshold,]), stats, plots))
}
