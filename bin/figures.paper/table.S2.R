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



#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")


# 2. groups' vector
groups <- c("blood", "iPS", "fibro_myoblasts", "muscle", "brain")

# 3. read df of group ELS, group-specific ELS, group-specific TSS-distal ELS

group.ELS <- c()
group.specific.ELS <- c()
group.specific.TSS.distal.ELS <- c()

for ( i in 1:5 ) {
  x <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".ELS.bed"),
                    stringsAsFactors = F, h=F, sep="\t")
  group.ELS <- c(group.ELS, nrow(x))
  
  y <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".specific.ELS.bed"),
                  stringsAsFactors = F, h=F, sep="\t")
  group.specific.ELS <- c(group.specific.ELS, nrow(y))
  
  z <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".specific.ELS.distal.2Kb.bed"),
                  stringsAsFactors = F, h=F, sep="\t")
  group.specific.TSS.distal.ELS <- c(group.specific.TSS.distal.ELS, nrow(z))
  
}

df <- data.frame(groups = groups,
                 group.ELS = group.ELS,
                 group.specific.ELS = group.specific.ELS,
                 group.specific.TSS.distal.ELS = group.specific.TSS.distal.ELS)

print(xtable(df), include.rownames=F)
