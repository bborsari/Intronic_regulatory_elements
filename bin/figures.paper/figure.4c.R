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
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo/")

# 2. read summary stats dataframe
# computed on TSS-distal ELS
m <- read.table("group.ELS/summary.stats.tsv", h=F, sep="\t")
colnames(m) <- c("total", "fraction", "type", "group")

# 3. read number of ELSs that are both intronic and exonic
x <- read.table("group.ELS/intronic_and_exonic.ELSs.tsv", h=F, sep="\t")
colnames(x) <- c("total", "fraction", "type", "group")


# 4. subtract from intronic and exonic those ELSs that are both intronic and exonic
m$fraction <- ifelse(m$type == "exonic", m$fraction - x$fraction, m$fraction)

m$percentage <- (m$fraction / m$total) * 100
m$group <- factor(m$group, levels = c("common", "stem_cells", "neural_progenitors", "differentiated_tissues"))
m$type <- gsub("intronic", "introns", m$type)

m$type <- gsub("exonic", "exons", m$type)
m$type <- gsub("intergenic", "intergenic regions", m$type)
m$type <- factor(m$type, levels = c("introns", "exons", "intergenic regions"))


# 3. make plot
palette = c("common" = "white",
            "stem_cells" = "gray",
            "differentiated_tissues" = "#41ab5d",
            "neural_progenitors" = "gold")



round_percent <- function(x) { 
  x <- x/sum(x)*100  # Standardize result
  res <- floor(x)    # Find integer bits
  rsum <- sum(res)   # Find out how much we are missing
  if(rsum<100) { 
    # Distribute points based on remainders and a random tie breaker
    o <- order(x%%1, sample(length(x)), decreasing=TRUE) 
    res[o[1:(100-rsum)]] <- res[o[1:(100-rsum)]]+1
  } 
  res 
}

round_percent(m[m$group=="differentiated_tissues", "fraction"] / m[m$group=="differentiated_tissues", "total"])



pdf("~/public_html/enhancers_neural_development/figures.paper/fig.4c.pdf", 
    width = 10, height = 4)
ggplot(m, aes(x=group, y=percentage, fill=group)) +
  geom_bar(stat = "identity", color = "black", alpha = .8) +
  geom_text(stat='identity', aes(label = round(percentage, 1)), vjust=-1,
            size = 5.5) +
  facet_grid(~type) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.background.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("proportion of ELSs (%)") +
  scale_fill_manual(values = palette) +
  ylim(0,100)
dev.off()
