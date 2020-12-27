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
library(cowplot)



#********
# BEGIN *
#********


# 1. number of introns per gene
source("/no_backup/rg/bborsari/projects/enhancers_neural_development/bin/figures.paper/figure.S3b1.R")
rm(m, w, x, groups, i)


# 2. median intron length per gene
source("/no_backup/rg/bborsari/projects/enhancers_neural_development/bin/figures.paper/figure.S3b2.R")
rm(m, w, x, groups, i)


# 3. density of ELSs per introns
source("/no_backup/rg/bborsari/projects/enhancers_neural_development/bin/figures.paper/figure.S3b3.R")
rm(m, x, groups, i)


# 4. merge plots
pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S3b.pdf",
    width = 14, height=4)
p <- plot_grid(plotlist = list(S3b1, S3b2, S3b3),
          nrow = 1, align = "h")
plot_grid(p, leg, nrow = 2, rel_heights = c(1, 0.1))
dev.off()