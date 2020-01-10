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
library(cowplot)


#************
# FUNCTIONS *
#************

my.function <- function(my.folder, my.table, my.xlab="") {
  
  # 1. read dataframes
  common.m <- read.table(paste0("intersection.with.adult/", my.folder, ".common.", my.table), h=F, sep="\t")
  colnames(common.m) <- c("region", "active_in_n_adult_tissues")
  common.m$group <- "common"
  
  differentiated_tissues.m <- read.table(paste0("intersection.with.adult/", my.folder, ".differentiated_tissues.", my.table), h=F, sep="\t")
  colnames(differentiated_tissues.m) <- c("region", "active_in_n_adult_tissues")
  differentiated_tissues.m$group <- "differentiated_tissues"
  
  neural_progenitors.m <- read.table(paste0("intersection.with.adult/", my.folder, ".neural_progenitors.", my.table), h=F, sep="\t")
  colnames(neural_progenitors.m) <- c("region", "active_in_n_adult_tissues")
  neural_progenitors.m$group <- "neural_progenitors"
  
  stem_cells.m <- read.table(paste0("intersection.with.adult/", my.folder, ".stem_cells.", my.table), h=F, sep="\t")
  colnames(stem_cells.m) <- c("region", "active_in_n_adult_tissues")
  stem_cells.m$group <- "stem_cells"
  
  
  # 2. merge datadframes
  merged.m <- rbind(common.m, differentiated_tissues.m, neural_progenitors.m, stem_cells.m)
  merged.m$group <- factor(merged.m$group, levels = c("stem_cells", "neural_progenitors", "differentiated_tissues", "common"))
  
  # 3. palette
  palette = c("differentiated_tissues" = "forestgreen",
              "neural_progenitors" = "goldenrod3",
              "mesoderm" = "#9e9ac8",
              "stem_cells" ="gray70",
              # "germ_layers" = "gray59",
              "common" = "white")
  
  p <- ggplot(merged.m, aes(x=active_in_n_adult_tissues, fill=group)) + 
    geom_density(color="black") +
    facet_wrap(~group, scales="free_y", nrow=1) +
    xlab(my.xlab) +
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = palette) +
    guides(fill=F)

  return(p)
  
}


#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo/")

lop <- list()

# 2. embryo: analysis all; adult: merged.table.tsv
lop[[1]] <- my.function(my.folder = "analysis.all", my.table = "merged.table.tsv")

# 3. embryo: analysis; adult: merged.table.subset.tsv
lop[[2]] <- my.function(my.folder = "analysis", my.table = "merged.table.subset.tsv", 
                        my.xlab = "# adult tissues in which the embryonic ELS are active")

# 4. plot
pdf("~/public_html/enhancers_neural_development/plots/fig.11.pdf", width = 13, height = 8)
plot_grid(plotlist = lop, nrow=2, ncol=1, align="v")
dev.off()
