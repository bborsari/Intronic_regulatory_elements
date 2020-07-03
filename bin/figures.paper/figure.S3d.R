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
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo")


# 2. groups' vector
groups <- c("stem_cells", "neural_progenitors", "differentiated_tissues", "common")


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

for ( i in 1:4 ) {
  
  m <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".genes.txt"),
                  stringsAsFactors = F, h=F, sep="\t")
  m <- merge(m, w, by.x = "V4", by.y = "gene_id")
  m$group <- groups[i]
  x <- rbind(x, m[, c("median", "group")])
  
}
x$group <- factor(x$group, levels = c("common",
                                      "stem_cells",
                                      "neural_progenitors",
                                      "differentiated_tissues"))


# 5. vector of % of intronic ELSs
y <- c(38.297872, 50.899743, 60.331633, 60.634648)
names(y) <- c("common",
              "stem_cells",
              "neural_progenitors",
              "differentiated_tissues")


# 6. test correlation between gene length and fraction of intronic ELSs
z <- x %>%
  group_by(group) %>%
  summarise(mean = mean(median), median = median(median), n = n())
stopifnot(identical(names(y), as.character(z$group)))

test <- cor.test(z$median, y, method = "pearson")
corr.v <- paste0("Pearson's R: ", 
                 round(test$estimate, 2), 
                 "; p = ",
                 format.pval(test$p.value, digits=2))


# 7. make plot
pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S3d.pdf",
    width = 4.5, height=3.5)
ggplot(x, aes(x=group, y=log10(median+1), fill=group)) +
  geom_violin(alpha=.65, aes(color = group)) +
  geom_boxplot(width=0.3, alpha = .85) +
  stat_compare_means(label.x = 1, label.y = 7.5) +
  annotate(geom="text", x=1.57, y=7, label=corr.v,
           color="black") +
  scale_fill_manual(values =c("common" = "white",
                              "stem_cells" = "gray",
                              "differentiated_tissues" = "#41ab5d",
                              "neural_progenitors" = "gold")) +
  scale_color_manual(values = c("common" = "#969696",
                                "stem_cells" = "white",
                                "differentiated_tissues" = "white",
                                "neural_progenitors" = "white")) +
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
  guides(fill=F, color=F) +
  ylab("log10(bp + 1)") +
  labs(title = "median intron length per gene") +
  ylim(0.5, 8)
dev.off()
