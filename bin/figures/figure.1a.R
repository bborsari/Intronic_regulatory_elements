.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(ggpubr)
library(cowplot)

setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/ELS.vs.psychENCODE/")



#************
# FUNCTIONS *
#************

map.pvals <- function(x){
  
  if (x < 0.001) {
    
    pval = "***"
    
  } else if (x < 0.01) {
    
    pval = "**"
    
  } else if (x < 0.05) {
    
    pval = "*"
     
  } else {
    
    pval = "ns"
    
  }
  
  return(pval)
  
}


#********
# BEGIN *
#********


df <- data.frame(intersecting = c(556460, 48917, 11046),
                 total = c(925190, 65983, 13073),
                 type = c("non-PFC", "also-PFC", "PFC-specific"))
df$type <- factor(df$type, levels = c("non-PFC", "also-PFC", "PFC-specific"))

# 1. compute p-values and OR of Fisher test

pvals <- c()
OR <- c()

## 1.1 non-PFC vs. also-PFC
m <- matrix(c(df[df$type=="non-PFC", "intersecting"],
              (df[df$type=="non-PFC", "total"] - df[df$type=="non-PFC", "intersecting"]),
              df[df$type=="also-PFC", "intersecting"],
              (df[df$type=="also-PFC", "total"] - df[df$type=="also-PFC", "intersecting"])),
            byrow = T, nrow = 2)

my.test <- fisher.test(m)
pvals <- c(pvals, my.test$p.value)
OR <- c(OR, my.test$estimate)


## 1.2 also-PFC vs. PFC-specific
m <- matrix(c(df[df$type=="also-PFC", "intersecting"],
              (df[df$type=="also-PFC", "total"] - df[df$type=="also-PFC", "intersecting"]),
              df[df$type=="PFC-specific", "intersecting"],
              (df[df$type=="PFC-specific", "total"] - df[df$type=="PFC-specific", "intersecting"])),
            byrow = T, nrow = 2)

my.test <- fisher.test(m)
pvals <- c(pvals, my.test$p.value)
OR <- c(OR, my.test$estimate)


## 1.3 non-PFC vs. PFC-specific
m <- matrix(c(df[df$type=="non-PFC", "intersecting"],
              (df[df$type=="non-PFC", "total"] - df[df$type=="non-PFC", "intersecting"]),
              df[df$type=="PFC-specific", "intersecting"],
              (df[df$type=="PFC-specific", "total"] - df[df$type=="PFC-specific", "intersecting"])),
            byrow = T, nrow = 2)

my.test <- fisher.test(m)
pvals <- c(pvals, my.test$p.value)
OR <- c(OR, my.test$estimate)


# 2. generate df with results from Fisher test
Fisher.res <- data.frame(pvals = pvals,
                         OR = OR,
                         combinations = c("non-PFC vs. also-PFC",
                                          "also-PFC vs. PFC-specific",
                                          "non-PFC vs. PFC-specific"))

Fisher.res$label <- sapply(Fisher.res$pvals, map.pvals)
Fisher.res$xmin <- c(1, 2, 1)
Fisher.res$xmax <- c(2, 3, 3)
Fisher.res$y_pos <- c(80, 90, 100)



# 3. make plot

pdf("~/public_html/enhancers_neural_development/plots/fig.1a.pdf", height = 6, width = 5)
ggplot(df, aes(x=type, y=(intersecting/total)*100, fill = type)) +
  geom_bar(stat="identity", width = .7, color = "black") +
  scale_fill_manual(values = c("#fed976", "#74c476", "#9e9ac8")) +
  guides(fill = F) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  ylim(0, 100) +
  ylab("% of REs intersecting introns") +
  geom_bracket(
    xmin = as.numeric(Fisher.res$xmin), 
    xmax = Fisher.res$xmax,
    y.position = Fisher.res$y_pos, 
    label = Fisher.res$label,
    tip.length = 0.05,
    inherit.aes = F,
    label.size = 7
  )

dev.off()
