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

#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/embryo/")


# 2. read dfs of embryo/adult intersection

df <- data.frame(stringsAsFactors = F)
for (x in c("common", "stem_cells", "neural_progenitors", "differentiated_tissues")) {
    
    tmp <- read.table(paste0("intersection.with.adult/", x,
                             ".specific.ELS.distal.2Kb.bed.tsv"), h=F, sep="\t")
    tmp$group <- x
    df <- rbind(df, tmp)
  
}


colnames(df)[1:2] <- c("ELS_id", "n_adult_samples")
df$group <- factor(df$group, levels = c("stem_cells", "neural_progenitors",
                                        "differentiated_tissues", "common"))

# 3. make plot
palette = c("common" = "white",
            "stem_cells" = "gray",
            "differentiated_tissues" = "#41ab5d",
            "neural_progenitors" = "gold")


pdf("~/public_html/enhancers_neural_development/figures.paper/fig.4d.pdf", 
    width = 9, height = 3.5)
ggplot(df, aes(x=n_adult_samples, fill = group)) +
  geom_histogram(color = "black", binwidth = 2) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  guides(fill = F) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        strip.background.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Number of active ELSs") +
  xlab("Number of adult samples")
dev.off()
