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
library(xtable)
library(ggpubr)
library(scales)




#********
# BEGIN *
#********


# 1. define working directory
setwd("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult")


# 2. read dfs of genic and intergenic ELSs
# and their degree of tissue-sharing
x <- read.table("ELS.sharing.genic.tsv", h=F, sep="\t")
y <- read.table("ELS.sharing.intergenic.tsv", h=F, sep="\t")


# 3. add info of genic/intergenic type

x$group <- "genic"
y$group <- "intergenic"


# 4. merge the 2 dfs
m <- rbind(x, y)
colnames(m)[1:5] <- c("chrom", "start", "end", "ELS", "n_samples")


# 5. keep only ELSs active in at least one tissue
m <- m[m$n_samples > 0,]


# 6. assign a bin to each ELS according to the number of samples
# it's active in
m$bins <- cut(m$n_samples, 
              breaks=c(seq(0, 60, by = 5)),
              labels = c("1-5", "6-10", "11-15", "16-20",
                         "21-25", "26-30", "31-35", "36-40",
                         "41-45", "46-50", "51-55", "56-60"))

# 7. summarize % of intergenic ELSs per nº of samples
p <- as.data.frame.matrix(table(m[, c("n_samples", "group")])) 
p <- as.data.frame(t(apply(p, 1, function(x){x <- x / sum(x)})))
p$n_samples <- as.numeric(rownames(p))


# 8. compute correlation (Pearson & Spearman)
# between % of intergenic ELSs and the nº of samples in which
# a ELS is active
test.p <- cor.test(p$intergenic, p$n_samples, method = "pearson")
test.s <- cor.test(p$intergenic, p$n_samples, method = "spearman")
corr.p <- paste0("Pearson's R: ", 
                 round(test.p$estimate, 2), 
                 "; p = ",
                 format.pval(test.p$p.value, digits=2))

corr.s <- paste0("Spearman's R: ", 
                 round(test.s$estimate, 2), 
                 "; p = ",
                 format.pval(test.s$p.value, digits=2))


# 9. make plot
pdf("~/public_html/enhancers_neural_development/figures.paper/fig.1a.pdf",
    width = 4, height=3.5)
ggplot(p, aes(x=n_samples, y=intergenic)) +
  geom_point(color = "black", size=3) +
  xlab("number of samples") +
  ylab("intergenic ELSs (%)") +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(labels = percent_format(), limits = c(0,1)) 
  # annotate(geom="text", x=25, y=0.95, label=corr.p,
  #          color="black") +
  # annotate(geom="text", x=27.5, y=0.85, label=corr.s,
  #          color="black")
dev.off()
  

