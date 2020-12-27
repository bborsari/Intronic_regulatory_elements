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
rownames(m) <- m$region
m$region <- NULL


# 2. read metadata
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
m <- m[, metadata$File_accession]
stopifnot(identical(colnames(m), metadata$File_accession))
colnames(m) <- metadata$Biosample_term_name



# 2. groups' vector
groups <- c("blood", "iPS", "fibro_myoblasts", "muscle", "brain")

# 3. read df of group ELS
x <- list()
for ( i in 1:5 ) {
  tmp <- read.table(paste0("group.ELS/", groups[i], "/", groups[i], ".ELS.bed"),
                       stringsAsFactors = F, h=F, sep="\t")
  x[[i]] <- paste(tmp$V1, tmp$V2, tmp$V3, tmp$V4, sep="_")
}

names(x) <- groups



# 3. groups2
groups2 <- list()
groups2$blood <- 1:5
groups2$iPS <- 6:7
groups2$fibro_myoblasts <- 8:10
groups2$muscle <- 11:18
groups2$brain <- 19:25


y <- data.frame(stringsAsFactors = F)

for ( i in 1:5 ) {
  
  for ( j in 1:5 ) {
    
    
    regions <- x[[i]]
    tmp <- m[regions, groups2[[j]]]
    print(nrow(tmp))
    tmp$sum <- apply(tmp, 1, sum)
    
    tmp2 <- data.frame(group1 = rep(groups[i], nrow(tmp)),
                       group2 = rep(groups[j], nrow(tmp)),
                       n = tmp$sum)
    
    y <- rbind(y, tmp2)
    
  }
  
}



y$group1 <- factor(y$group1, levels = c("blood", "iPS",
                                        "fibro_myoblasts",
                                        "muscle", "brain"))
y$group2 <- factor(y$group2, levels = c("blood", "iPS",
                                        "fibro_myoblasts",
                                        "muscle", "brain"))


palette = c("brain" = "#EEEE00",
            "blood" = "#FF00BB",
            "muscle" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue")

y <- y[y$group1 != y$group2, ]

pdf("~/public_html/enhancers_neural_development/figures.paper/fig.S1b.pdf",
    height = 7.5, width = 8)
ggplot(y[y$n > 0, ], aes(x=n, fill=group2)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_grid(group1~., scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 8, by = 1),
                   labels = as.character(seq(0, 8, 1))) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title = element_text(size =20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_text(size=20, angle = 0),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size=20)) +
  xlab("Number of samples") +
  ylab("Number of active ELSs") +
  labs(title="Sharing of tissue-active ELSs among clusters")
dev.off()





