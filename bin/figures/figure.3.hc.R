.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

# run as an interactive job

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
library(dendextend)


#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/")

# 2. import binary dataframe
m <- as.data.frame(fread("merged.table.tsv"))
rownames(m) <- m$region
m$region <- NULL

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
filtered.m <- m[, colnames(m) %in% (metadata[metadata$File_accession %in% selected.adult.samples, "File_accession"])]
colnames(filtered.m) <- metadata[metadata$File_accession %in% selected.adult.samples, "Biosample_term_name"]


# 7. define palette
palette = c("brain" = "gold",
            "blood" = "#FF61CC",
            "mesoderm" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue")


# 8. assign emvryological groups
groups <- c()
filtered.metadata <- metadata[metadata$File_accession %in% selected.adult.samples, ]

for (i in 1:nrow(filtered.metadata)) {
  
  if (filtered.metadata[i, "Biosample_term_name"] %in% c("B_cell",
                                                         "CD14-positive_monocyte",
                                                         "T-cell",
                                                         "natural_killer_cell",
                                                         "peripheral_blood_mononuclear_cell")) {
    
    groups <- c(groups, "blood") 
    
  } else if (filtered.metadata[i, "Biosample_term_name"] %in% c("gastrocnemius_medialis",
                                                                "skeletal_muscle_tissue",
                                                                "rectal_smooth_muscle_tissue",
                                                                "right_cardiac_atrium",
                                                                "stomach_smooth_muscle",
                                                                "subcutaneous_abdominal_adipose_tissue",
                                                                "muscle_layer_of_duodenum",
                                                                "vagina")) {
    groups <- c(groups, "mesoderm")
    
  } else if (filtered.metadata[i, "Biosample_term_name"] %in% c("angular_gyrus",
                                                                "caudate_nucleus",
                                                                "cingulate_gyrus",
                                                                "layer_of_hippocampus",
                                                                "middle_frontal_area_46",
                                                                "substantia_nigra",
                                                                "temporal_lobe")) {
    
    groups <- c(groups, "brain")
    
  } else if (filtered.metadata[i, "Biosample_term_name"] %in% c("myotube",
                                                                "fibroblast_of_lung",
                                                                "skeletal_muscle_myoblast")) {
    
    groups <- c(groups, "fibro_myoblasts")
    
    
  } else {
    
    groups <- c(groups, "iPS")
    
  }
  
}

filtered.metadata$groups <- groups
rownames(filtered.metadata) <- filtered.metadata$Biosample_term_name


my.ann_colors <- list(groups = palette)

out <- pheatmap(filtered.m, 
         show_rownames = F,
         cluster_rows = F, 
         cluster_cols = T, 
         clustering_distance_cols = "binary",
         clustering_method = "ward.D",
         annotation_col = filtered.metadata[, "groups", drop=F],
         annotation_colors = my.ann_colors)
dev.off()



pdf("~/public_html/enhancers_neural_development/plots/fig.3.hc.pdf")
plot(out$tree_col)
dev.off()
