.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(data.table)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
library(ggplot2)
library(ggrepel)

#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/")

# 2. read dataframe
m <- read.table("analysis/intersection.intronic.tsv", h=F, sep="\t")
colnames(m) <- c("total", "fraction", "group", "distance")

m$percentage <- (m$fraction / m$total) * 100
m$group <- factor(m$group, levels = c("common", "iPS", "fibro_myoblasts",
                                      "blood", "mesoderm", "brain"))
# 3. make plot
palette = c("brain" = "gold",
            "endoderm" = "brown",
            "blood" = "#FF61CC",
            "mesoderm" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue",
            "common" = "white")


pdf("~/public_html/enhancers_neural_development/plots/fig.4.pdf", width = 10, height = 4.5)
ggplot(m, aes(x=group, y=percentage, fill=group)) +
  geom_bar(stat = "identity", color = "black", alpha = .8) +
  facet_grid(~distance) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        legend.position = "bottom",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("% of ELS overlapping introns") +
  scale_fill_manual(values = palette) +
  ylim(0,100)
dev.off()
