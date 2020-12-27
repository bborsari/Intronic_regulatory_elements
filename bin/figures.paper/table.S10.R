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
library(scales)

#********
# BEGIN *
#********

# 1. source previous script to compute % 
# of intronic / exonic / intergenic ELSs
source("/no_backup/rg/bborsari/projects/enhancers_neural_development/bin/figures.paper/figure.4c.R")


# 2. perform pairwise Fisher test
# between tissue-specific and common ELSs
pvals <- c()
OR <- c()
CI <- c()
groups <- c()
types <- c()

for (type in c("introns", "exons", "intergenic regions")) {
  
  tmp <- m[m$type == type, ]
  tmp$total <- tmp$total - tmp$fraction
  
  for (group in (c("stem_cells", "neural_progenitors", "differentiated_tissues"))) {
    
    tmp2 <- tmp[tmp$group %in% c(group, "common"), 1:2]
    ftest <- fisher.test(tmp2)
    pvals <- c(pvals, ftest$p.value)
    OR <- c(OR, ftest$estimate)
    CI <- c(CI, paste(round(ftest$conf.int[[1]], 2), round(ftest$conf.int[[2]], 2), sep ="-"))
    groups <- c(groups, group)
    types <- c(types, type)
    
  }
  
  
}


# 3. build dataframe for table
y <- data.frame(pvals = pvals,
                groups = groups,
                OR = OR,
                CI = CI,
                types = types)

# 4. compute FDR values
y$FDR <- p.adjust(y$pvals, method = "BH")
y$FDR <- scientific(y$FDR, digits = 2)

# 5. make table
y <- y[, c("types", "groups", "FDR", "OR", "CI")]
colnames(y)[1:2] <- c("genomic location", "group")
y$OR <- round(y$OR, 2)

print(xtable(y), include.rownames=F)
