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
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo/")

# 2. import binary dataframe
m <- as.data.frame(fread("merged.table.tsv"))
rownames(m) <- m$region
m$region <- NULL

# 5. read metadata
metadata <- read.table("5-group.ccREs.hg19.bigBed.txt", h=T, sep="\t")
selected.embryo.samples <- c("ENCFF840ANN",
                             "ENCFF903RGX",
                             "ENCFF543DVJ",
                             "ENCFF800YES",
                             "ENCFF941JIE",
                             "ENCFF198WHL",
                             "ENCFF093BQM",
                             "ENCFF059PHA",
                             "ENCFF292NZP",
                             "ENCFF281QON",
                             "ENCFF477EUQ",
                             "ENCFF112ZGF",
                             "ENCFF455CQW",
                             "ENCFF138OGZ",
                             "ENCFF593TNG",
                             "ENCFF376XBS",
                             "ENCFF051OUV",
                             "ENCFF021HBJ",
                             "ENCFF505OUS",
                             "ENCFF086FKD",
                             "ENCFF205SDB",
                             "ENCFF180QLH",
                             "ENCFF332EYK",
                             "ENCFF138DOQ",
                             "ENCFF620BVM",
                             "ENCFF250CGY")

# 6. filter for manually selected samples
filtered.m <- m[, colnames(m) %in% (metadata[metadata$File_accession %in% selected.embryo.samples, "File_accession"])]
colnames(filtered.m) <- metadata[metadata$File_accession %in% selected.embryo.samples, "Biosample_term_name"]


# 7. define palette
palette = c("differentiated_tissues" = "forestgreen",
            "neural_progenitors" = "goldenrod3",
            "mesoderm" = "#9e9ac8",
            "stem_cells" ="gray70",
            "germ_layers" = "gray59",
            "common" = "white")


# 8. assign embryological groups
groups <- c()
filtered.metadata <- metadata[metadata$File_accession %in% selected.embryo.samples, ]

for (i in 1:nrow(filtered.metadata)) {
  
  if (filtered.metadata[i, "Biosample_term_name"] %in% c("adrenal_gland",
                                                         "large_intestine",
                                                         "small_intestine",
                                                         "muscle_of_trunk",
                                                         "muscle_of_leg",
                                                         "stomach",
                                                         "hepatocyte",
                                                         "thymus",
                                                         "fibroblast_of_lung",
                                                         "smooth_muscle_cell")) {
    
    groups <- c(groups, "differentiated_tissues") 
    
  } else if (filtered.metadata[i, "Biosample_term_name"] %in%  c("neural_cell",
                                                                 "neural_progenitor_cell",
                                                                 "neural_stem_progenitor_cell",
                                                                 "neuroepithelial_stem_cell",
                                                                 "radial_glial_cell",
                                                                 "mid-neurogenesis_radial_glial_cells")) {
    groups <- c(groups, "neural_progenitors")
    
  } else if (filtered.metadata[i, "Biosample_term_name"] %in% c("H1",
                                                                "H9",
                                                                "HUES48",
                                                                "HUES6",
                                                                "HUES64")) {
    
    groups <- c(groups, "stem_cells")
    
    
  } else {
    
    groups <- c(groups, "germ_layers")
    
  }
  
}

filtered.metadata$groups <- groups
filtered.metadata$Biosample_term_name <- as.character(filtered.metadata$Biosample_term_name)
filtered.metadata[filtered.metadata$File_accession == "ENCFF505OUS", "Biosample_term_name"] <- "H9_1"
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



pdf("~/public_html/enhancers_neural_development/plots/fig.7.hc.pdf")
plot(out$tree_col)
dev.off()
