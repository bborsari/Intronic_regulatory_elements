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
library(plotly)
library(scatterplot3d)
library(VennDiagram)
library(xtable)
library(ggpubr)



#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")


# 2. groups' vector
groups <- c("blood", "iPS", "fibro_myoblasts", "muscle", "brain", "common")


# 3. read df with intron length per gene
# and compute mean/median intron length per gene
w <- read.table("gencode.v19.intron.sizes.tsv", h=F, sep="\t")


colnames(w) <- c("gene_id", "intron_length")
w <- w %>%
  group_by(gene_id) %>%
  summarise(mean = mean(intron_length), median = median(intron_length), n = n())


# 4. read dfs of genes intersecting with group-specific ccREs
# and compute number of introns per gene
x <- data.frame(stringsAsFactors = F)

for ( i in 1:6 ) {
  
  m <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".genes.txt"),
                  stringsAsFactors = F, h=F, sep="\t")
  
  m <- merge(m, w, by.x = "V4", by.y = "gene_id")
  m$group <- groups[i]
  x <- rbind(x, m[, c("median", "group")])
  
}
x$group <- factor(x$group, levels = c("common",
                                      "iPS",
                                      "fibro_myoblasts",
                                      "blood",
                                      "muscle",
                                      "brain"))


# 5. make plot
# pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S2a2.pdf",
#     width = 4.5, height=3.5)
S2a2 <- ggplot(x, aes(x=group, y=log10(median+1), fill=group)) +
  geom_violin(alpha=.65, aes(color = group)) +
  geom_boxplot(width=0.3, alpha = .85) +
  stat_compare_means(label.x = 1.5, label.y = 7.5, size =6) +
  scale_fill_manual(values =c("brain" = "gold",
                              "blood" = "#FF61CC",
                              "muscle" = "#9e9ac8",
                              "iPS" ="gray",
                              "fibro_myoblasts" = "lightblue",
                              "common" = "white")) +
  scale_color_manual(values = c("brain" = "white",
                                "blood" = "white",
                                "muscle" = "white",
                                "iPS" ="white",
                                "fibro_myoblasts" = "white",
                                "common" = "#969696")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        plot.title = element_text(size=20, hjust=.5),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=F, color=F) +
  ylab("bp (log10)") +
  labs(title = "median intron length") +
  ylim(0.5, 8)
# dev.off()
