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
source("/no_backup/rg/bborsari/projects/enhancers_neural_development/bin/figures.paper/figure.S1a.R")

# 5. selected adult samples
selected.adult.samples <- c("ENCFF379TAE",
                            "ENCFF967MJU",
                            "ENCFF098NHL",
                            "ENCFF529UWB",
                            "ENCFF509DPX",
                            "ENCFF863OGG",
                            "ENCFF311MNY",
                            "ENCFF093MDL",
                            "ENCFF278RUJ",
                            "ENCFF726JTT",
                            "ENCFF725QLM",
                            "ENCFF862BGI",
                            "ENCFF904XYE",
                            "ENCFF942KAC",
                            "ENCFF508GKP",
                            "ENCFF494WCN",
                            "ENCFF159NZA",
                            "ENCFF070EXF",
                            "ENCFF233VRB",
                            "ENCFF810IQU",
                            "ENCFF120MMC",
                            "ENCFF495RTY",
                            "ENCFF037UZZ",
                            "ENCFF920QRH",
                            "ENCFF231KWX")

# 6. subset fit
fit2 <- fit[fit$File_accession %in% selected.adult.samples, ]


# 7. assign germ-layer groups
groups <- c()
for (i in 1:nrow(fit2)) {
  
  if (fit2[i, "tissues"] %in% c("B_cell",
                               "CD14-positive_monocyte",
                               "T-cell",
                               "natural_killer_cell",
                               "peripheral_blood_mononuclear_cell")) {
    
    groups <- c(groups, "blood") 
    
  } else if (fit2[i, "tissues"] %in% c("gastrocnemius_medialis",
                                      "skeletal_muscle_tissue",
                                      "rectal_smooth_muscle_tissue",
                                      "right_cardiac_atrium",
                                      "stomach_smooth_muscle",
                                      "subcutaneous_abdominal_adipose_tissue",
                                      "muscle_layer_of_duodenum",
                                      "vagina")) {
    groups <- c(groups, "mesoderm")
    
  } else if (fit2[i, "tissues"] %in% c("angular_gyrus",
                                      "caudate_nucleus",
                                      "cingulate_gyrus",
                                      "layer_of_hippocampus",
                                      "middle_frontal_area_46",
                                      "substantia_nigra",
                                      "temporal_lobe")) {
    
    groups <- c(groups, "brain")
    
  } else if (fit2[i, "tissues"] %in% c("myotube",
                                      "fibroblast_of_lung",
                                      "skeletal_muscle_myoblast")) {
    
    groups <- c(groups, "fibro_myoblasts")
    
    
  } else {
    
    groups <- c(groups, "iPS")
    
  }
  
}

fit2$groups <- groups

# 8. define color palette for germ-layer groups
palette = c("brain" = "#EEEE00",
            "blood" = "#FF00BB",
            "mesoderm" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue",
            "common" = "white")


# 9. make plot
pdf("~/public_html/enhancers_neural_development/figures.paper/fig.1b.pdf",
    height = 5, width = 6, useDingbats = F)
s3d <- scatterplot3d(fit2[, c(2,4,3)], pch = 16, 
                     color = palette[fit2$groups],
                     col.axis = "lightgray", 
                     grid = F,
                     ylim = c(-0.3, 0.3))
text(s3d$xyz.convert(fit2[, c(2,4,3)]), labels = fit2$File_accession,
     cex= 0.7)
dev.off()


# plot_ly(fit2, x = ~V1, y = ~V3, z = ~V2, 
#         color = ~groups,
#         text = ~tissues,
#         colors = palette)
