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

# 2. read summary stats dataframe
# computed on TSS-distal ELS
m <- read.table("group.ELS/summary.stats.tsv", h=F, sep="\t")
colnames(m) <- c("total", "fraction", "type", "group")

m$percentage <- (m$fraction / m$total) * 100
m$group <- factor(m$group, levels = c("common", "iPS", "fibro_myoblasts",
                                      "blood", "muscle", "brain"))
m$type <- gsub("intronic", "introns", m$type)
m$type <- gsub("exonic", "exons", m$type)
m$type <- gsub("intergenic", "intergenic regions", m$type)
m$type <- factor(m$type, levels = c("introns", "exons", "intergenic regions"))


# 3. make plot
palette = c("brain" = "#EEEE00",
            "blood" = "#FF00BB",
            "muscle" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue",
            "common" = "white")



pdf("~/public_html/enhancers_neural_development/figures.paper/fig.2a.pdf", 
    width = 9, height = 4)
ggplot(m, aes(x=group, y=percentage, fill=group)) +
  geom_bar(stat = "identity", color = "black", alpha = .8) +
  geom_text(stat='identity', aes(label = paste0(round(percentage, 1), "%")), vjust=-1) +
  facet_grid(~type) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        strip.background.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("% of ELSs") +
  scale_fill_manual(values = palette) +
  ylim(0,100)
dev.off()
