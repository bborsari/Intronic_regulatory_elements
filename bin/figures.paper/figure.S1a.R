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


#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")

# 2. import binary dataframe
# m <- as.data.frame(fread("merged.table.tsv"))
# rownames(m) <- m$region
# m$region <- NULL

# 3. compute distance matrix
# t.m <- t(m)
# dist.t.m <- dist(t.m, method = "binary")
# dist.t.m.mat <- as.matrix(dist.t.m)

# 4. load pre-computed distance matrix
load("dist.Rdata")

# 5. read metadata
metadata <- read.table("5-group.ccREs.hg19.bigBed.txt", h=T, sep="\t")

# 6. retrieve distance object
distance.object <- as.dist(dist.t.m.mat)

# 7. perform mds
fit <- cmdscale(distance.object, eig=TRUE, k=3) # k is the number of dim
stopifnot(identical(colnames(dist.t.m.mat), rownames(fit$points))) # view results
fit <- as.data.frame(fit$points)
fit$tissues <- metadata[metadata$File_accession %in% rownames(dist.t.m.mat), "Biosample_term_name"]


# 8. retrieve file_id
fit$File_accession <- rownames(fit)


# 9. read color palette
palette <- read.delim("palette.txt", h=T, sep="\t", stringsAsFactors = F)


# 10. add to fit info of colors
fit <- merge(fit, palette[, c("File_accession", "GTEx.color")], 
             by = "File_accession")


# 11. make 3D plot
pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S1a.pdf",
    height = 5, width = 6, useDingbats = F)
s3d <- scatterplot3d(fit[, c(2,4,3)], pch = 16, 
                     color = fit[, "GTEx.color"], 
                     col.axis = "lightgray", 
                     grid = F,
                     ylim = c(-0.3, 0.3))
text(s3d$xyz.convert(fit[, c(2,4,3)]), labels = fit$File_accession,
     cex= 0.3)
dev.off()

