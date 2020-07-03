setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/")


#************
# LIBRARIES *
#************

library(pheatmap)
library(grid)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)



## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_75 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 75, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_75 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 75, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_75",
                  ns=asNamespace("pheatmap"))


common.H3K27ac <- read.table("sandra.H3K27ac.txt", h=F, sep="\t")
colnames(common.H3K27ac) <- "enhancer"
common.H3K27ac$status <- "shared_H3K27ac"

common.H3K27ac.HiC <- read.table("sandra.H3K27ac.HiC.txt", h=F, sep="\t")
colnames(common.H3K27ac.HiC) <- "enhancer"
common.H3K27ac.HiC$status2 <- "shared_H3K27ac_HiC"


###-------
### TFs
###-------

# read table 
merged.TFs <- read.table("TFs.chromatin.remodellers/analysis/all.files.tsv", h=F, sep="\t",
                         stringsAsFactors = F)
merged.TFs$enhancer <- paste(merged.TFs$V1, merged.TFs$V2, merged.TFs$V3, sep="_")

# retrieve set of experiments
experiments.TFs <- unique(merged.TFs$V15)
experiments.TFs.df <- data.frame(experiments = experiments.TFs)
experiments.TFs.df <- 
  experiments.TFs.df %>% separate(experiments, c("target", "biosample", "file_id"), ";")

# write.table(experiments.TFs.df, "~/public_html/enhancers_neural_development/TFs/table.experiments.tsv",
#             quote=F, row.names = F, col.names = T, sep="\t")

# pdf("~bborsari/public_html/enhancers_neural_development/plots/table.ChIP-seq.TFs.pdf", width=8)
# grid.table(experiments.TFs.df, rows = NULL)
# dev.off()




###------
# retrieve presence/absence of peaks
###------

# dataframe where to store presence/absence of peaks
final.TFs <- data.frame(enhancer=unique(merged.TFs$enhancer), stringsAsFactors = F)
final.TFs <- final.TFs[order(final.TFs$enhancer), ,drop=F]

for (exp in experiments.TFs) {
  
  tmp <- merged.TFs[merged.TFs$V15 == exp, ]
  tmp <- tmp[order(tmp$enhancer, tmp$V14), ]
  tmp <- tmp[match(unique(tmp[, "enhancer"]), tmp[, "enhancer"]),]
  tmp$V14 <- ifelse((tmp$V14>0), 1, 0)
  final.TFs <- merge(final.TFs, tmp[, c("V14", "enhancer")], by="enhancer")
  colnames(final.TFs)[ncol(final.TFs)] <- exp
  
}

final.TFs <- merge(final.TFs, common.H3K27ac, by="enhancer", all = T)
final.TFs$status <- ifelse(is.na(final.TFs$status), "not_shared_H3K27ac", "shared_H3K27ac")
stopifnot(sum(final.TFs$status == "shared_H3K27ac") == 393)

final.TFs <- merge(final.TFs, common.H3K27ac.HiC, by="enhancer", all = T)
final.TFs$status2 <- ifelse(is.na(final.TFs$status2), "not_shared_H3K27ac_HiC", "shared_H3K27ac_HiC")
stopifnot(sum(final.TFs$status2 == "shared_H3K27ac_HiC") == 123)

rownames(final.TFs) <- final.TFs$enhancer
final.TFs$enhancer <- NULL


# pdf("~bborsari/public_html/enhancers_neural_development/heatmap.ChIP-seq.TFs.pdf", width=14)
pheatmap(as.matrix(final.TFs[, 1:23]), clustering_method = "ward.D", 
         clustering_distance_rows = "binary",
         clustering_distance_cols = "binary",
         legend = F,
         show_rownames = F,
         fontsize = 14, 
         annotation_row = final.TFs[, 24:25])
# dev.off()


###------
# retrieve only regions with at least one binding event across experiments
###------

final.TFs$my_sum <- apply(final.TFs[, 1:23], 1, sum)
final.TFs.sub <- final.TFs[final.TFs$my_sum >0, 1:25]

pdf("~/public_html/enhancers_neural_development/plots/heatmap.ChIP-seq.TFs.at.least.one.binding.pdf", width=14)
p.TFs <- pheatmap(as.matrix(final.TFs.sub[, 1:23]), clustering_method = "ward.D", 
         clustering_distance_rows = "binary",
         clustering_distance_cols = "binary",
         legend = F,
         show_rownames = F,
         fontsize = 14, 
         annotation_row = final.TFs.sub[, 24:25])
dev.off()

###------
# check elbow plot to decide optimal number of clusters (8)
###------

## 1. compute distance
X.TFs <- dist(final.TFs.sub[, 1:23], method="binary")

## 2. retrieve wssval
# Auxiliary function to compute WSS
wss <- function(d) {
  sum(scale(d, scale = FALSE)^2)
}
wssval.TFs <- c()
for(i in 1:20){
  # Compute hierarchical clustering
  hc <- hclust(X.TFs, method = "ward.D")
  
  # Cut tree
  clusts <- cutree(hc, k = i)
  
  # Split data per clusters and calculate WSS for each cluster, and sum
  spl <- split(X.TFs, clusts)
  wssval.TFs[i] <- sum(sapply(spl, wss))
}


# 3. Generate plot
# Create plot data
plot.TFs.elbow <- data.frame(wss = wssval.TFs,
                             nb_clust = seq_along(wssval.TFs))
# Create differences between clusters to show on the plot
lags <- round((wssval.TFs - lag(wssval.TFs, 1))[-1], digits = 3)
# Plot
pdf("~/public_html/enhancers_neural_development/plots/elbow.plot.TFs.with.CTCF.pdf")
ggplot(plot.TFs.elbow, aes(x = nb_clust, y = wss)) + 
  geom_point(size = 3) +
  geom_line() + 
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5)) +
  scale_x_continuous(breaks = seq(1, 20, 1)) +
  xlab("Number of Clusters") +
  ylab("Within groups sum of squares") +
  annotate("text", x = seq(1, 19, 1) + 0.6, y = wssval.TFs[-20] + 0.1, label = lags) +
  labs(title = "TFs - elbow plot")
dev.off()


###----
# cut df with appropriate number of clusters
###----
cluster_row_TFs = as.data.frame( cutree( p.TFs$tree_row, k = 8 ) )
colnames(cluster_row_TFs) <- "hierarchical_clustering"
cluster_row_TFs$region <- rownames(cluster_row_TFs)
cluster_row_TFs <- cluster_row_TFs %>% separate(region, c("chrom", "start", "end"), "_")
cluster_row_TFs$hc <- cluster_row_TFs$hierarchical_clustering
cluster_row_TFs$hierarchical_clustering <- NULL
write.table(cluster_row_TFs, "~/public_html/enhancers_neural_development/TFs/hc.tree.with.CTCF.tsv",
            col.names = T, row.names = F, sep="\t", quote=F)




###-----
### repeat elbow plot removing experiments on CTCF
###-----

pdf("~/public_html/enhancers_neural_development/plots/heatmap.ChIP-seq.TFs.at.least.one.binding.pdf.no.CTCF.pdf", width=14)
p.TFs.no.CTCF <- pheatmap(as.matrix(final.TFs.sub[, c(3:6, 8:22)]), clustering_method = "ward.D", 
                          clustering_distance_rows = "binary",
                          clustering_distance_cols = "binary",
                          legend = F,
                          show_rownames = F,
                          fontsize = 14, 
                          annotation_row = final.TFs.sub[, 24:25])
dev.off()

###------
# check elbow plot to decide optimal number of clusters (8)
###------

## 1. compute distance
X.TFs.no.CTCF <- dist(final.TFs.sub[, c(3:6, 8:22)], method="binary")

## 2. retrieve wssval
# Auxiliary function to compute WSS
wss <- function(d) {
  sum(scale(d, scale = FALSE)^2)
}
wssval.TFs.no.CTCF <- c()
for(i in 1:20){
  # Compute hierarchical clustering
  hc <- hclust(X.TFs.no.CTCF, method = "ward.D")
  
  # Cut tree
  clusts <- cutree(hc, k = i)
  
  # Split data per clusters and calculate WSS for each cluster, and sum
  spl <- split(X.TFs.no.CTCF, clusts)
  wssval.TFs.no.CTCF[i] <- sum(sapply(spl, wss))
}


# 3. Generate plot
# Create plot data
plot.TFs.no.CTCF.elbow <- data.frame(wss = wssval.TFs.no.CTCF,
                                     nb_clust = seq_along(wssval.TFs.no.CTCF))
# Create differences between clusters to show on the plot
lags <- round((wssval.TFs.no.CTCF - lag(wssval.TFs.no.CTCF, 1))[-1], digits = 3)
# Plot
pdf("~/public_html/enhancers_neural_development/plots/elbow.plot.TFs.no.CTCF.pdf")
ggplot(plot.TFs.no.CTCF.elbow, aes(x = nb_clust, y = wss)) + 
  geom_point(size = 3) +
  geom_line() + 
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5)) +
  scale_x_continuous(breaks = seq(1, 20, 1)) +
  xlab("Number of Clusters") +
  ylab("Within groups sum of squares") +
  annotate("text", x = seq(1, 19, 1) + 0.6, y = wssval.TFs.no.CTCF[-20] + 0.1, label = lags) +
  labs(title = "TFs (no CTCF) - elbow plot")
dev.off()


###----
# cut df with appropriate number of clusters
###----
cluster_row_TFs.no.CTCF = as.data.frame( cutree( p.TFs.no.CTCF$tree_row, k = 8 ) )
colnames(cluster_row_TFs.no.CTCF) <- "hierarchical_clustering"
cluster_row_TFs.no.CTCF$region <- rownames(cluster_row_TFs.no.CTCF)
cluster_row_TFs.no.CTCF <- cluster_row_TFs.no.CTCF %>% separate(region, c("chrom", "start", "end"), "_")
cluster_row_TFs.no.CTCF$hc <- cluster_row_TFs.no.CTCF$hierarchical_clustering
cluster_row_TFs.no.CTCF$hierarchical_clustering <- NULL
write.table(cluster_row_TFs.no.CTCF, "~/public_html/enhancers_neural_development/TFs/hc.tree.no.CTCF.tsv",
            col.names = T, row.names = F, sep="\t", quote=F)










###-------
### histone marks (replicated peaks)
###-------

# read table 
merged.histone.marks.replicated <- 
  read.table("histone.marks/replicated_peaks/analysis/all.files.tsv", h=F, sep="\t",
                         stringsAsFactors = F)
merged.histone.marks.replicated$enhancer <- paste(merged.histone.marks.replicated$V1, 
                                                  merged.histone.marks.replicated$V2, 
                                                  merged.histone.marks.replicated$V3, 
                                                  sep="_")

# retrieve set of experiments
experiments.histone.marks.replicated <- sort(unique(merged.histone.marks.replicated$V15))
experiments.histone.marks.replicated.df <- data.frame(experiments = experiments.histone.marks.replicated)
experiments.histone.marks.replicated.df <- 
  experiments.histone.marks.replicated.df %>% separate(experiments, c("target", "biosample", "file_id"), ";")

# write.table(experiments.histone.marks.replicated.df, "~/public_html/enhancers_neural_development/histone.marks/replicated_peaks/table.experiments.tsv",
#             quote=F, row.names = F, col.names = T, sep="\t")


# pdf("~bborsari/public_html/enhancers_neural_development/plots/table.ChIP-seq.histone.marks.replicated.pdf", width=8, height=10)
# grid.table(experiments.histone.marks.replicated.df, rows = NULL)
# dev.off()




# dataframe where to store presence/absence of peaks
final.histone.marks.replicated <- 
  data.frame(enhancer=unique(merged.histone.marks.replicated$enhancer), 
             stringsAsFactors = F)
final.histone.marks.replicated <- 
  final.histone.marks.replicated[order(final.histone.marks.replicated$enhancer), ,drop=F]

# retrieve presence/absence of peaks
for (exp in experiments.histone.marks.replicated) {
  
  tmp <- merged.histone.marks.replicated[merged.histone.marks.replicated$V15 == exp, ]
  tmp <- tmp[order(tmp$enhancer, tmp$V14), ]
  tmp <- tmp[match(unique(tmp[, "enhancer"]), tmp[, "enhancer"]),]
  tmp$V14 <- ifelse((tmp$V14>0), 1, 0)
  final.histone.marks.replicated <- merge(final.histone.marks.replicated, tmp[, c("V14", "enhancer")], by="enhancer")
  colnames(final.histone.marks.replicated)[ncol(final.histone.marks.replicated)] <- exp
  
}

final.histone.marks.replicated <- merge(final.histone.marks.replicated, common.H3K27ac, by="enhancer", all = T)
final.histone.marks.replicated$status <- ifelse(is.na(final.histone.marks.replicated$status), "not_shared_H3K27ac", "shared_H3K27ac")
stopifnot(sum(final.histone.marks.replicated$status == "shared_H3K27ac") == 393)

final.histone.marks.replicated <- merge(final.histone.marks.replicated, common.H3K27ac.HiC, by="enhancer", all = T)
final.histone.marks.replicated$status2 <- ifelse(is.na(final.histone.marks.replicated$status2), "not_shared_H3K27ac_HiC", "shared_H3K27ac_HiC")
stopifnot(sum(final.histone.marks.replicated$status2 == "shared_H3K27ac_HiC") == 123)

rownames(final.histone.marks.replicated) <- final.histone.marks.replicated$enhancer
final.histone.marks.replicated$enhancer <- NULL


pdf("~bborsari/public_html/enhancers_neural_development/heatmap.ChIP-seq.histone.marks.replicated.peaks.pdf", width=14)
pheatmap(as.matrix(final.histone.marks.replicated[, 1:28]), clustering_method = "ward.D", 
         cluster_cols = F,
         clustering_distance_rows = "binary",
         legend = F,
         show_rownames = F,
         annotation_row = final.histone.marks.replicated[, 29:30])
dev.off()






###-------
### histone marks (stable peaks)
###-------

# read table 
merged.histone.marks.stable <- 
  read.table("histone.marks/stable_peaks/analysis/all.files.tsv", h=F, sep="\t",
             stringsAsFactors = F)
merged.histone.marks.stable$enhancer <- paste(merged.histone.marks.stable$V1, 
                                                  merged.histone.marks.stable$V2, 
                                                  merged.histone.marks.stable$V3, 
                                                  sep="_")

# retrieve set of experiments
experiments.histone.marks.stable <- sort(unique(merged.histone.marks.stable$V15))
experiments.histone.marks.stable.df <- data.frame(experiments = experiments.histone.marks.stable)
experiments.histone.marks.stable.df <- 
  experiments.histone.marks.stable.df %>% separate(experiments, c("target", "biosample", "file_id"), ";")

# write.table(experiments.histone.marks.stable.df, "~/public_html/enhancers_neural_development/histone.marks/stable_peaks/table.experiments.tsv",
#             quote=F, row.names = F, col.names = T, sep="\t")


# pdf("~bborsari/public_html/enhancers_neural_development/table.ChIP-seq.histone.marks.stable.pdf", width=8, height=10)
# grid.table(experiments.histone.marks.stable.df, rows = NULL)
# dev.off()




# dataframe where to store presence/absence of peaks
final.histone.marks.stable <- 
  data.frame(enhancer=unique(merged.histone.marks.stable$enhancer), 
             stringsAsFactors = F)
final.histone.marks.stable <- 
  final.histone.marks.stable[order(final.histone.marks.stable$enhancer), ,drop=F]

# retrieve presence/absence of peaks
for (exp in experiments.histone.marks.stable) {
  
  tmp <- merged.histone.marks.stable[merged.histone.marks.stable$V15 == exp, ]
  tmp <- tmp[order(tmp$enhancer, tmp$V14), ]
  tmp <- tmp[match(unique(tmp[, "enhancer"]), tmp[, "enhancer"]),]
  tmp$V14 <- ifelse((tmp$V14>0), 1, 0)
  final.histone.marks.stable <- merge(final.histone.marks.stable, tmp[, c("V14", "enhancer")], by="enhancer")
  colnames(final.histone.marks.stable)[ncol(final.histone.marks.stable)] <- exp
  
}

final.histone.marks.stable <- merge(final.histone.marks.stable, common.H3K27ac, by="enhancer", all = T)
final.histone.marks.stable$status <- ifelse(is.na(final.histone.marks.stable$status), "not_shared_H3K27ac", "shared_H3K27ac")
stopifnot(sum(final.histone.marks.stable$status == "shared_H3K27ac") == 393)

final.histone.marks.stable <- merge(final.histone.marks.stable, common.H3K27ac.HiC, by="enhancer", all = T)
final.histone.marks.stable$status2 <- ifelse(is.na(final.histone.marks.stable$status2), "not_shared_H3K27ac_HiC", "shared_H3K27ac_HiC")
stopifnot(sum(final.histone.marks.stable$status2 == "shared_H3K27ac_HiC") == 123)

rownames(final.histone.marks.stable) <- final.histone.marks.stable$enhancer
final.histone.marks.stable$enhancer <- NULL


pdf("~bborsari/public_html/enhancers_neural_development/heatmap.ChIP-seq.histone.marks.stable.peaks.pdf", width=14)
pheatmap(as.matrix(final.histone.marks.stable[, 1:100]), clustering_method = "ward.D", 
         cluster_cols = F,
         clustering_distance_rows = "binary",
         legend = F,
         show_rownames = F,
         annotation_row = final.histone.marks.stable[, 101:102],
         annotation_legend = F)
dev.off()
