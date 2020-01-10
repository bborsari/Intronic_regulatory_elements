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


# 6. filter for manually selected samples
filtered.dist.t.m.mat <- dist.t.m.mat[rownames(dist.t.m.mat) %in% (metadata[metadata$File_accession %in% selected.adult.samples, 
                                                                            "File_accession"]), 
                                      colnames(dist.t.m.mat) %in% (metadata[metadata$File_accession %in% selected.adult.samples, 
                                                                            "File_accession"])]

# 7. retrieve distance object
distance.object <- as.dist(filtered.dist.t.m.mat)

# 8. perform mds
fit <- cmdscale(distance.object, eig=TRUE, k=3) # k is the number of dim
stopifnot(identical(colnames(filtered.dist.t.m.mat), rownames(fit$points))) # view results
fit <- as.data.frame(fit$points)
fit$tissues <- metadata[metadata$File_accession %in% rownames(filtered.dist.t.m.mat), "Biosample_term_name"]

# 9. assign emvryological groups
groups <- c()
for (i in 1:nrow(fit)) {
  
  if (fit[i, "tissues"] %in% c("B_cell",
                               "CD14-positive_monocyte",
                               "T-cell",
                               "natural_killer_cell",
                               "peripheral_blood_mononuclear_cell")) {
    
    groups <- c(groups, "blood") 
    
  } else if (fit[i, "tissues"] %in% c("gastrocnemius_medialis",
                                      "skeletal_muscle_tissue",
                                      "rectal_smooth_muscle_tissue",
                                      "right_cardiac_atrium",
                                      "stomach_smooth_muscle",
                                      "subcutaneous_abdominal_adipose_tissue",
                                      "muscle_layer_of_duodenum",
                                      "vagina")) {
    groups <- c(groups, "mesoderm")
    
  } else if (fit[i, "tissues"] %in% c("angular_gyrus",
                                      "caudate_nucleus",
                                      "cingulate_gyrus",
                                      "layer_of_hippocampus",
                                      "middle_frontal_area_46",
                                      "substantia_nigra",
                                      "temporal_lobe")) {
    
    groups <- c(groups, "brain")
    
  } else if (fit[i, "tissues"] %in% c("myotube",
                                      "fibroblast_of_lung",
                                      "skeletal_muscle_myoblast")) {
    
    groups <- c(groups, "fibro_myoblasts")
    
    
  } else {
    
    groups <- c(groups, "iPS")
    
  }
  
}

fit$groups <- groups

# 10. define palette
palette = c("brain" = "gold",
            "endoderm" = "brown",
            "blood" = "#FF61CC",
            "mesoderm" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue",
            "common" = "white")


# 11. plot
mds.plot.adult <- plot_ly(fit, x = ~V1, y = ~V2, z = ~V3, 
                          color = ~groups,
                          text = ~tissues,
                          colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(mds.plot.adult), "~/public_html/enhancers_neural_development/plots/fig.3.mds.html")

