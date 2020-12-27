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
library(scales)


#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo")

# 2. read binary table
m <- fread("merged.table.subset.tsv")
m <- as.data.frame(m)
rownames(m) <- m$region
m$region <- NULL


# 3. read metadata
metadata <- read.table("5-group.ccREs.hg19.bigBed.txt", h=T, sep="\t",
                       stringsAsFactors = F)
metadata <- metadata[metadata$File_accession %in% colnames(m), ]
stopifnot(identical(colnames(m), metadata$File_accession))
metadata$Biosample_term_name <- factor(metadata$Biosample_term_name,
                                       levels = c("HUES6",
                                                  "HUES64",
                                                  "HUES48",
                                                  "H9",
                                                  "H9.1",
                                                  "H1",
                                                  "neuroepithelial_stem_cell",
                                                  "radial_glial_cell",
                                                  "neural_progenitor_cell",
                                                  "mid-neurogenesis_radial_glial_cells",
                                                  "neural_stem_progenitor_cell",
                                                  "neural_cell",
                                                  "smooth_muscle_cell",
                                                  "thymus",
                                                  "adrenal_gland",
                                                  "fibroblast_of_lung",
                                                  "muscle_of_trunk",
                                                  "muscle_of_leg",
                                                  "stomach",
                                                  "hepatocyte",
                                                  "large_intestine",
                                                  "small_intestine"))

metadata <- metadata[order(metadata$Biosample_term_name), ]
metadata$group <- c(rep("stem_cells", 6), 
                    rep("neural_progenitors", 6),
                    rep("differentiated_tissues", 10))

m <- m[, metadata$File_accession]
stopifnot(identical(colnames(m), metadata$File_accession))
colnames(m) <- metadata$Biosample_term_name

metadata[metadata$File_accession == "ENCFF505OUS", "Biosample_term_name"] <- "H9.1"


# 4. groups' vector
groups <- c("stem_cells", "neural_progenitors", "differentiated_tissues")

# 5. read df of group-specific ELS
x <- list()
for ( i in 1:3 ) {
  tmp <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".specific.ELS.distal.2Kb.bed"),
                    stringsAsFactors = F, h=F, sep="\t")
  x[[i]] <- paste(tmp$V1, tmp$V2, tmp$V3, tmp$V4, sep="_")
}

names(x) <- groups



# 6. define groups2
groups2 <- list()
groups2$stem_cells <- 1:6
groups2$neural_progenitors <- 7:12
groups2$differentiated_tissues <- 13:22



# 7. compute degree of contamination
# as the non-group tissues in which the ccREs are active
df <- data.frame(stringsAsFactors = F)

for ( i in 1:3 ) {
  
  tmp <- m[x[[i]], -c(groups2[[i]])]

  y <- c()
  for ( j in 1:nrow(tmp) ){
    
    if (sum(tmp[j, ]) == 1) {
      
      z <- colnames(tmp)[which(tmp[j, ] == 1)]
      y <- c(y, metadata[metadata$Biosample_term_name == z, "group"])
      
    } else {
      
      y <- c(y, "within-cluster")
      
    }
    
    
  }
  
  y <- data.frame(table(y))
  y$group <- groups[i]
  
  df <- rbind(df, y)
  
}




# 8. define color palette
palette = c("stem_cells" = "gray",
            "differentiated_tissues" = "#41ab5d",
            "neural_progenitors" = "gold",
            "within-cluster" = "#525252")


# 9. make plot
df$group <- factor(df$group, levels = c("stem_cells", "differentiated_tissues", "neural_progenitors"))
df$y <- factor(df$y, levels = c("stem_cells", "neural_progenitors", "differentiated_tissues", "within-cluster"))


pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S3a.pdf",
    height = 5, width = 8)
ggplot(df, aes(x=group, y=Freq, fill=y)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title = element_text(size =20),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 30, vjust = .6),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_text(size=20, angle = 0),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  scale_y_continuous(labels = percent_format()) +
  ylab("proportion (%)") +
  xlab("groups")+
  scale_x_discrete(labels = c("stem\ncells", "differentiated\ntissues", "neural\nprogenitors"))

dev.off()
