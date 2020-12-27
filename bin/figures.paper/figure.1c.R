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


#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")


# 2. read binary table
m <- fread("merged.table.subset.tsv")
m <- as.data.frame(m)


# 3. compute number of samples in which the ccRE is active
m$sum <- apply(m[, 2:ncol(m)], 1, sum)


# 4. keep only ccREs active in at least two tissues
m <- m[m$sum > 1, ]
rownames(m) <- m$region
m$region <- NULL
m$sum <- NULL


# 5. read metadata
metadata <- read.table("5-group.ccREs.hg19.bigBed.txt", h=T, sep="\t",
                       stringsAsFactors = F)
metadata <- metadata[metadata$File_accession %in% colnames(m), ]
stopifnot(identical(colnames(m), metadata$File_accession))


# 6. reorder metadata according to table S1 order
t1.order <- c("B_cell",
              "T-cell",
              "natural_killer_cell",
              "peripheral_blood_mononuclear_cell",
              "CD14-positive_monocyte",
              "iPS-18a",
              "iPS-20b",
              "fibroblast_of_lung",
              "skeletal_muscle_myoblast",
              "myotube",
              "muscle_layer_of_duodenum",
              "subcutaneous_abdominal_adipose_tissue",
              "rectal_smooth_muscle_tissue",
              "vagina",
              "stomach_smooth_muscle",
              "skeletal_muscle_tissue",
              "right_cardiac_atrium",
              "gastrocnemius_medialis",
              "middle_frontal_area_46",
              "caudate_nucleus",
              "angular_gyrus",
              "layer_of_hippocampus",
              "substantia_nigra",
              "temporal_lobe",
              "cingulate_gyrus")
metadata$Biosample_term_name <- factor(metadata$Biosample_term_name,
                                       levels = t1.order)

metadata <- metadata[order(metadata$Biosample_term_name), ]
metadata$cluster <- c(rep("blood", 5), rep("iPS", 2),
                    rep("fibro_myoblasts", 3), 
                    rep("muscle", 8),
                    rep("brain", 7))
rownames(metadata) <- metadata$Biosample_term_name

# 7. reorder binary table according to metadata file
m <- m[, metadata$File_accession]
stopifnot(identical(colnames(m), metadata$File_accession))


# 8. convert colnames of binary table from file_ids 
# to biosample term names
colnames(m) <- metadata$Biosample_term_name


# 9. compute minimum overlap of active ccREs between 
# two samples
x <- data.frame(stringsAsFactors = F)

for ( i in 1:25 ) {
  
  for ( j in 1:25 ) {
    
    tmp <- melt(table(m[, c(i, j)]))
    colnames(tmp) <- c("from", "to", "n")
    tmp$ex1 <- colnames(m)[i]
    tmp$ex2 <- colnames(m)[j]
    n_ex1 <- sum(tmp[tmp$from == 1, "n"])
    tmp$n <- tmp$n / n_ex1
    
    x <- rbind(x, tmp[tmp$from == 1 & tmp$to ==1, ])
    
  }
  
}

x <- x[, c("ex1", "ex2", "n")]


# 10. convert x into matrix
y <- dcast(x, formula = ex1~ex2)
rownames(y) <- y$ex1
y$ex1 <- NULL


# 11. reorder matrix for heatmap
# z <- c("B_cell",
#        "T-cell",
#        "natural_killer_cell",
#        "peripheral_blood_mononuclear_cell",
#        "CD14-positive_monocyte",
#        "iPS-18a",
#        "iPS-20b",
#        "myotube",
#        "skeletal_muscle_myoblast",
#        "fibroblast_of_lung",
#        "muscle_layer_of_duodenum",
#        "subcutaneous_abdominal_adipose_tissue",
#        "rectal_smooth_muscle_tissue",
#        "vagina",
#        "stomach_smooth_muscle",
#        "skeletal_muscle_tissue",
#        "right_cardiac_atrium",
#        "gastrocnemius_medialis",
#        "middle_frontal_area_46",
#        "caudate_nucleus",
#        "angular_gyrus",
#        "layer_of_hippocampus",
#        "substantia_nigra",
#        "temporal_lobe",
#        "cingulate_gyrus")


z <- c("iPS-18a",
       "iPS-20b",
       "peripheral_blood_mononuclear_cell",
       "T-cell",
       "natural_killer_cell",
       "B_cell",
       "CD14-positive_monocyte",
       "middle_frontal_area_46",
       "angular_gyrus",
       "temporal_lobe",
       "caudate_nucleus",
       "cingulate_gyrus",
       "layer_of_hippocampus",
       "substantia_nigra",
       "rectal_smooth_muscle_tissue",
       "stomach_smooth_muscle",
       "muscle_layer_of_duodenum",
       "vagina",
       "subcutaneous_abdominal_adipose_tissue",
       "right_cardiac_atrium",
       "skeletal_muscle_tissue",
       "gastrocnemius_medialis",
       "myotube",
       "skeletal_muscle_myoblast",
       "fibroblast_of_lung")


y <- y[z, z]

labs <- c("50", "52", "11", "9", "10", "5", "14",
          "54", "56", "59", "55", "60",
          "57", "58", "39", "41", "36",
          "40", "38", "43", "42", "44", "24",
          "23", "22")

# 11. heatmap
pdf("~/public_html/enhancers_neural_development/figures.paper/fig.1c.pdf",
    height = 10, width = 12)
pheatmap(y*100, cluster_cols = T, cluster_rows = F,
         clustering_method = "complete",
         border_color = NA,
         cellwidth = 20,
         cellheight = 20,
         fontsize = 20,
         labels_row = labs,
         labels_col = labs,
         annotation_row = metadata[z, "cluster", drop=F],
         annotation_col = metadata[z, "cluster", drop=F],
         annotation_colors = list(cluster = c("brain" = "#EEEE00",
                                            "blood" = "#FF00BB",
                                            "muscle" = "#9e9ac8",
                                            "iPS" ="gray",
                                            "fibro_myoblasts" = "lightblue")))
dev.off()



# old order:
# c("B_cell",
#   "T-cell",
#   "natural_killer_cell",
#   "peripheral_blood_mononuclear_cell",
#   "CD14-positive_monocyte",
#   "iPS-18a",
#   "iPS-20b",
#   "myotube",
#   "skeletal_muscle_myoblast",
#   "fibroblast_of_lung",
#   "subcutaneous_abdominal_adipose_tissue",
#   "muscle_layer_of_duodenum",
#   "gastrocnemius_medialis",
#   "skeletal_muscle_tissue",
#   "right_cardiac_atrium",
#   "rectal_smooth_muscle_tissue",
#   "stomach_smooth_muscle",
#   "vagina",
#   "angular_gyrus",
#   "cingulate_gyrus",
#   "temporal_lobe",
#   "middle_frontal_area_46",
#   "layer_of_hippocampus",
#   "caudate_nucleus",
#   "substantia_nigra")