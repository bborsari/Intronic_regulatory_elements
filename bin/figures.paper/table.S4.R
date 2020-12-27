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
groups <- c("blood", "iPS", "fibro_myoblasts", "muscle", "brain", "common")


# 3. list of genes with 
# intronic or exonic ccREs
x.intronic <- c()
x.exonic <- c()
x.intronic.exonic <- c()

for ( i in 1:6 ) {
  
  tmp.intronic <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".genes.intronic.txt"),
                             stringsAsFactors = F, h=F, sep="\t")
  tmp.intronic <- tmp.intronic$V1
  
  tmp.exonic <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".genes.exonic.txt"),
                             stringsAsFactors = F, h=F, sep="\t")
  tmp.exonic <- tmp.exonic$V1
  
  x.intronic <- c(x.intronic,
                  length(setdiff(tmp.intronic, tmp.exonic)))
  
  x.exonic <- c(x.exonic,
                length(setdiff(tmp.exonic, tmp.intronic)))
  
  x.intronic.exonic <- c(x.intronic.exonic,
                         length(intersect(tmp.exonic, tmp.intronic)))
  
  
}


# 4. prepare summary table
df <- data.frame(groups = groups,
                 g_intronic_ELS = x.intronic,
                 g_exonic_ELS = x.exonic,
                 g_intronic_exonic_ELS = x.intronic.exonic)
df$total <- apply(df[, 2:4], 1, sum)


df$g_intronic_ELS <- paste0(df$g_intronic_ELS, " (", 
                            round(((df$g_intronic_ELS/df$total)*100), 2), "%)")
df$g_exonic_ELS <- paste0(df$g_exonic_ELS, " (", 
                          round(((df$g_exonic_ELS/df$total)*100), 2), "%)")
df$g_intronic_exonic_ELS <- paste0(df$g_intronic_exonic_ELS, " (", 
                                   round(((df$g_intronic_exonic_ELS/df$total)*100), 2), "%)")


# 5. print table in latex format
print(xtable(df), include.rownames=F)
