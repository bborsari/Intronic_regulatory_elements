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
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")

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
                                       levels = c("B_cell",
                                                  "T-cell",
                                                  "natural_killer_cell",
                                                  "peripheral_blood_mononuclear_cell",
                                                  "CD14-positive_monocyte",
                                                  "iPS-18a",
                                                  "iPS-20b",
                                                  "myotube",
                                                  "skeletal_muscle_myoblast",
                                                  "fibroblast_of_lung",
                                                  "subcutaneous_abdominal_adipose_tissue",
                                                  "muscle_layer_of_duodenum",
                                                  "gastrocnemius_medialis",
                                                  "skeletal_muscle_tissue",
                                                  "right_cardiac_atrium",
                                                  "rectal_smooth_muscle_tissue",
                                                  "stomach_smooth_muscle",
                                                  "vagina",
                                                  "angular_gyrus",
                                                  "cingulate_gyrus",
                                                  "temporal_lobe",
                                                  "middle_frontal_area_46",
                                                  "layer_of_hippocampus",
                                                  "caudate_nucleus",
                                                  "substantia_nigra"))

metadata <- metadata[order(metadata$Biosample_term_name), ]
metadata$group <- c(rep("blood", 5), rep("iPS", 2),
                    rep("fibro_myoblasts", 3), 
                    rep("muscle", 8),
                    rep("brain", 7))

m <- m[, metadata$File_accession]
stopifnot(identical(colnames(m), metadata$File_accession))
colnames(m) <- metadata$Biosample_term_name



# 4. groups' vector
groups <- c("blood", "iPS", "fibro_myoblasts", "muscle", "brain")

# 5. read df of group-specific ELS
x <- list()
for ( i in 1:5 ) {
  tmp <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".specific.ELS.distal.2Kb.bed"),
                    stringsAsFactors = F, h=F, sep="\t")
  x[[i]] <- paste(tmp$V1, tmp$V2, tmp$V3, tmp$V4, sep="_")
}

names(x) <- groups



# 6. define groups2
groups2 <- list()
groups2$blood <- 1:5
groups2$iPS <- 6:7
groups2$fibro_myoblasts <- 8:10
groups2$muscle <- 11:18
groups2$brain <- 19:25


# 7. compute degree of contamination
# as the non-group tissues in which the ccREs are active
df <- data.frame(stringsAsFactors = F)

for ( i in 1:5 ) {
  
  tmp <- m[x[[i]], -c(groups2[[i]])]
  
  y <- c()
  for ( j in 1:nrow(tmp) ){
    
    if (sum(tmp[j, ]) == 1) {
      
      z <- colnames(tmp)[which(tmp[j, ] == 1)]
      y <- c(y, metadata[metadata$Biosample_term_name == z, "group"])
      
    } else {
      
      y <- c(y, "within-group")
      
    }
    
    
  }
  
  y <- data.frame(table(y))
  y$group <- groups[i]
  
  df <- rbind(df, y)
  
}




# 8. define color palette
palette = c("brain" = "#EEEE00",
            "blood" = "#FF00BB",
            "muscle" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue",
            "within-group" = "#525252")


# 9. make plot
df$group <- factor(df$group, levels = c("iPS", "fibro_myoblasts", "brain", "blood", "muscle"))
df$y <- factor(df$y, levels = c("iPS", "blood", "brain", "fibro_myoblasts", "muscle", "within-group"))


pdf("~/public_html/enhancers_neural_development/figures.paper/fig.1d.pdf",
    height = 4, width = 5)
ggplot(df, aes(x=group, y=Freq, fill=y)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title = element_text(size =12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 30, vjust = .5),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_text(size=13, angle = 0),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  scale_y_continuous(labels = percent_format()) +
  ylab("frequency") +
  xlab("groups")

dev.off()
