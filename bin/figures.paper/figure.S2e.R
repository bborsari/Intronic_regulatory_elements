# .libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


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
# setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")


# 2. groups' vector
groups <- c("blood", "iPS", "fibro_myoblasts", "muscle", "brain", "common")


# 3. read dfs of introns intersecting with group-specific ccREs
x <- data.frame(stringsAsFactors = F)

for ( i in 1:6 ) {
  
  # m <- read.table(paste0("group.ELS/", 
  #                        groups[i], "/", 
  #                        groups[i],
  #                        ".specific.ELS.distal.2Kb.per.intron.bed"),
  #                 stringsAsFactors = F, h=F, sep="\t")
  
  m <- read.table(url(paste0("https://public-docs.crg.es/rguigo/Data/bborsari/enhancers_neural_development/ccREs/hg19/adult/group.ELS/", 
                             groups[i], "/", groups[i], ".specific.ELS.distal.2Kb.per.intron.bed")),
                  stringsAsFactors = F, h=F, sep="\t")
  
  
  m$group <- groups[i]
  x <- rbind(x, m[, c(8,9)])
  
}

colnames(x)[1] <- "n_ELS"

x$group <- factor(x$group, levels = c("common",
                                      "iPS",
                                      "fibro_myoblasts",
                                      "blood",
                                      "muscle",
                                      "brain"))


# 4. vector of % of intronic ELSs
y <- c(35.57692, 54.93742, 55.03123, 56.92008, 71.62068, 73.99605)
names(y) <- c("common",
              "iPS",
              "fibro_myoblasts",
              "blood",
              "muscle",
              "brain")


# 5. test correlation between gene length and fraction of intronic ELSs
z <- x %>%
  group_by(group) %>%
  summarise(mean = mean(n_ELS), median = median(n_ELS), n = n())
stopifnot(identical(names(y), as.character(z$group)))

test <- cor.test(z$median, y, method = "pearson")
corr.v <- paste0("Pearson's R: ", 
                 round(test$estimate, 2), 
                 "; p = ",
                 format.pval(test$p.value, digits=2))


# 7. make plot
# pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S2e.pdf",
#    width = 4.5, height=3.5)
ggplot(x, aes(x=group, y=n_ELS, fill=group)) +
  geom_violin(alpha=.65, aes(color = group), scale = "width") +
  geom_boxplot(width=0.3, alpha = .7) +
  stat_compare_means(label.x = 1.5, label.y = 2) +
  annotate(geom="text", x=2.4, y=60, label=corr.v,
           color="black") +
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
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        plot.title = element_text(size=14, hjust=.5),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(color=F, fill=F) +
  ylab("n") +
  labs(title = "number of ELSs per intron") +
  scale_y_continuous(trans = "log10")
# dev.off()
