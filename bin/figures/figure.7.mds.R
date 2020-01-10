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
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo")

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
                             "ENCFF180QLH")
                             # "ENCFF332EYK",
                             # "ENCFF138DOQ",
                             # "ENCFF620BVM",
                             # "ENCFF250CGY")


# 6. filter for manually selected samples
filtered.dist.t.m.mat <- dist.t.m.mat[rownames(dist.t.m.mat) %in% (metadata[metadata$File_accession %in% selected.embryo.samples, 
                                                                            "File_accession"]), 
                                      colnames(dist.t.m.mat) %in% (metadata[metadata$File_accession %in% selected.embryo.samples, 
                                                                            "File_accession"])]

# 7. retrieve distance object
distance.object <- as.dist(filtered.dist.t.m.mat)

# 8. perform mds
fit <- cmdscale(distance.object, eig=TRUE, k=3) # k is the number of dimensions
stopifnot(identical(colnames(filtered.dist.t.m.mat), rownames(fit$points))) # view results
fit <- as.data.frame(fit$points)
fit$tissues <- metadata[metadata$File_accession %in% rownames(filtered.dist.t.m.mat), "Biosample_term_name"]

# 9. assign emvryological groups
groups <- c()
for (i in 1:nrow(fit)) {
  
  if (fit[i, "tissues"] %in% c("adrenal_gland",
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
    
  } else if (fit[i, "tissues"] %in% c("neural_cell",
                                      "neural_progenitor_cell",
                                      "neural_stem_progenitor_cell",
                                      "neuroepithelial_stem_cell",
                                      "radial_glial_cell",
                                      "mid-neurogenesis_radial_glial_cells")) {
    groups <- c(groups, "neural_progenitors")
    
  } else if (fit[i, "tissues"] %in% c("H1",
                                      "H9",
                                      "HUES48",
                                      "HUES6",
                                      "HUES64")) {
    
    groups <- c(groups, "stem_cells")
    
    
  } #else {
    
    #groups <- c(groups, "germ_layers")
    
  #}
  
}

fit$groups <- groups

# 10. define palette
palette = c("differentiated_tissues" = "forestgreen",
            "neural_progenitors" = "goldenrod3",
            "mesoderm" = "#9e9ac8",
            "stem_cells" ="gray70",
            # "germ_layers" = "gray59",
            "common" = "white")


# 11. plot
mds.plot.embryo <- plot_ly(fit, x = ~V1, y = ~V2, z = ~V3, 
                          color = ~groups,
                          text = ~tissues,
                          colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(mds.plot.embryo), "~/public_html/enhancers_neural_development/plots/fig.7.mds.without.germ.layers.html")

