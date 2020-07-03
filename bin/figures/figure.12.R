#************
# LIBRARIES *
#************

library(ggplot2)

palette = c("brain" = "gold",
            "endoderm" = "brown",
            "blood" = "#FF61CC",
            "mesoderm" = "#9e9ac8",
            "iPS" ="gray",
            "fibro_myoblasts" = "lightblue",
            "common" = "white")


#************
# FUNCTIONS *
#************

# function 1 - gene length of genes intersecting with ccREs
my.function1 <- function(folder) {
  
  
  my.list <- list()
  
  # 1. import bedtools intersect table of genes intersecting group-specific ccREs
  for (i in 1:6) {
    
    my.list[[i]] <- read.table(paste0("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/", 
                                      folder, 
                                      "/",
                                      groups[i],
                                      "/stats/genes.intersecting.",
                                      groups[i],
                                      ".ccREs.tsv"),
                               h=F,
                               sep="\t")
    
  }
  
  
  # 2. for each group, retrieve length of genes intersecting with the ccREs
  my.df <- data.frame(stringsAsFactors = F)
  
  for ( i in 1:6 ) {
    
    tmp <- my.list[[i]]
    tmp2 <- data.frame(gene_length=tmp[, 3] - tmp[, 2],
                       group = rep(groups[i], nrow(tmp)))
    my.df <- rbind(my.df, tmp2)
    
  }
  
  my.df$group <- factor(my.df$group, levels = c("common", "iPS", "fibro_myoblasts",
                                                "blood","mesoderm", "brain"))
  
  p <- ggplot(my.df, aes(x=group, y=log10(gene_length+1), fill=group)) + 
    geom_boxplot() +
    stat_compare_means(comparisons = list(c("common", "iPS"),
                                          c("common", "fibro_myoblasts"),
                                          c("common", "blood"),
                                          c("common", "mesoderm"),
                                          c("common", "brain")),
                       method = "wilcox.test",
                       size = 4.5) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    ylab("gene length - log10 (bp+1)") +
    scale_fill_manual(values=palette) +
    guides(fill=F)
  
  return(p)
  
}


# function 2 - number of introns per gene intersecting group-specific ccREs
my.function2 <- function(folder) {
  
  my.list <- list()
  
  # 1. import bedtools intersect table of genes intersecting group-specific ccREs
  for (i in 1:6) {
    
    my.list[[i]] <- read.table(paste0("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/", 
                                      folder, 
                                      "/",
                                      groups[i],
                                      "/stats/genes.intersecting.",
                                      groups[i],
                                      ".ccREs.tsv"),
                               h=F,
                               sep="\t")
    
  }
  
  
  # 2. import table with number of introns per gene
  m <- read.table("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/gencode.v19.number.introns.per.gene.tsv",
  h=F, sep="\t")
  
  
  # 3. for each group, retrieve number of introns per each gene intersecting ccREs
  my.df <- data.frame(stringsAsFactors = F)
  
  for ( i in 1:6 ) {
    
    
    tmp <- my.list[[i]]
    tmp2 <- merge(tmp, m, by.x = "V4", by.y = "V1")
    tmp2$group <- groups[i]
    my.df <- rbind(my.df, tmp2[, c(ncol(tmp2)-1, ncol(tmp2))])
    
  }
  
  colnames(my.df)[1] <- "number_introns"
  
  my.df$group <- factor(my.df$group, levels = c("common", "iPS", "fibro_myoblasts",
                                                "blood","mesoderm", "brain"))
  
  p <- ggplot(my.df, aes(x=group, y=number_introns, fill=group)) + geom_boxplot() +
    stat_compare_means(comparisons = list(c("common", "iPS"),
                                          c("common", "fibro_myoblasts"),
                                          c("common", "blood"),
                                          c("common", "mesoderm"),
                                          c("common", "brain")),
                       method = "wilcox.test",
                       size = 4.5) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    ylab("n. of introns per gene") +
    scale_fill_manual(values=palette) +
    guides(fill=F)
  

  return(p)
  
}


# function 3 - length of all introns from genes intersecting group-specific ccREs
my.function3 <- function(folder) {
  
  
  my.list <- list()
  
  # 1. import bedtools intersect table of genes intersecting group-specific ccREs
  for (i in 1:6) {
    
    my.list[[i]] <- read.table(paste0("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/", 
                                      folder, 
                                      "/",
                                      groups[i],
                                      "/stats/genes.intersecting.",
                                      groups[i],
                                      ".ccREs.tsv"),
                               h=F,
                               sep="\t")
    
  }
 
  # 2. import table with length of non-redundant introns per gene
  m <- read.table("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/gencode.v19.intron.sizes.tsv",
                   h=F, sep="\t")
  
  my.df <- data.frame(stringsAsFactors = F)
  
  for ( i in 1:6 ) {
    
    tmp <- my.list[[i]]
    tmp2 <- m[m$V1 %in% tmp$V4, ]
    tmp2$group <- groups[i]
    my.df <- rbind(my.df, tmp2[, 2:3])
    
  }
  
  colnames(my.df)[1] <- "intron_length"
  
  my.df$group <- factor(my.df$group, levels = c("common", "iPS", "fibro_myoblasts",
                                                "blood","mesoderm", "brain"))
  
  p <- ggplot(my.df, aes(x=group, y=log10(intron_length+1), fill=group)) + geom_boxplot() +
    stat_compare_means(comparisons = list(c("common", "iPS"),
                                          c("common", "fibro_myoblasts"),
                                          c("common", "blood"),
                                          c("common", "mesoderm"),
                                          c("common", "brain")),
                       method = "wilcox.test",
                       size = 4.5)  +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    ylab("intron length - log10 (bp+1)") +
    scale_fill_manual(values=palette) +
    guides(fill=F) 
  
  
  return(p)

}

# function 4 - number of ccREs intersecting introns
my.function4 <- function(folder) {
  
  my.list <- list()
  
  # 1. import bedtools intersect table of genes intersecting group-specific ccREs
  for (i in 1:6) {
    
    my.list[[i]] <- read.table(paste0("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/", 
                                      folder, 
                                      "/",
                                      groups[i],
                                      "/stats/",
                                      groups[i],                                      
                                      ".ccREs.per.intron.tsv"),
                               h=F,
                               sep="\t")
    
  }
  
  my.df <- data.frame(stringsAsFactors = F)
  
  for ( i in 1:6 ) {
    
    tmp <- my.list[[i]]
    tmp$group <- groups[i]
    my.df <- rbind(my.df, tmp[, c(ncol(tmp)-1, ncol(tmp))])
    
  }
  
  colnames(my.df)[1] <- "number_enhancers"
  
  my.df$group <- factor(my.df$group, levels = c("common", "iPS", "fibro_myoblasts",
                                                "blood","mesoderm", "brain"))
  
  p <- ggplot(my.df, aes(x=group, y=log10(number_enhancers+1), fill=group)) + geom_boxplot() +
    stat_compare_means(comparisons = list(c("common", "iPS"),
                                          c("common", "fibro_myoblasts"),
                                          c("common", "blood"),
                                          c("common", "mesoderm"),
                                          c("common", "brain")),
                       method = "wilcox.test",
                       size = 4.5)  +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    ylab("n. of ccREs per intron - log10(n+1)") +
    scale_fill_manual(values=palette) +
    guides(fill=F) 
  
  return(p)
  
  
}



# function 5 - length of introns intersecting ccREs
my.function5 <- function(folder) {
  
  my.list <- list()
  
  # 1. import bedtools intersect table of genes intersecting group-specific ccREs
  for (i in 1:6) {
    
    my.list[[i]] <- read.table(paste0("/no_backup/rg/bborsari/projects/enhancers_neural_development/5-group.ccREs/hg19/adult/", 
                                      folder, 
                                      "/",
                                      groups[i],
                                      "/stats/",
                                      groups[i],                                      
                                      ".ccREs.per.intron.tsv"),
                               h=F,
                               sep="\t")
    
  }
  
  my.df <- data.frame(stringsAsFactors = F)
  
  for ( i in 1:6 ) {
    
    tmp <- my.list[[i]]
    tmp$group <- groups[i]
    tmp$length <- tmp[, 3] - tmp[, 2]
    my.df <- rbind(my.df, tmp[, c("group", "length")])
    
  }
  
  my.df$group <- factor(my.df$group, levels = c("common", "iPS", "fibro_myoblasts",
                                                "blood","mesoderm", "brain"))
  

  p <- ggplot(my.df, aes(x=group, y=log10(length+1), fill=group)) + geom_boxplot() +
    stat_compare_means(comparisons = list(c("common", "iPS"),
                                          c("common", "fibro_myoblasts"),
                                          c("common", "blood"),
                                          c("common", "mesoderm"),
                                          c("common", "brain")),
                       method = "wilcox.test",
                       size = 4.5)  +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    ylab("intersecting intron length - log10 (bp+1)") +
    scale_fill_manual(values=palette) +
    guides(fill=F) 
  
  return(p)
  
  
}


#********
# BEGIN *
#********

groups <- c("blood", "brain", "mesoderm", "fibro_myoblasts", "iPS", "common")
lop.analysis <- list()
lop.analysis.all <- list()


# 1. gene length of genes intersecting with ccREs
lop.analysis[[1]] <- my.function1(folder = "analysis")
lop.analysis.all[[1]] <- my.function1(folder = "analysis.all")

# 2. number of introns per gene intersecting with ccREs
lop.analysis[[2]] <- my.function2(folder="analysis")
lop.analysis.all[[2]] <- my.function2(folder="analysis.all")

# 3. length of all introns from genes intersecting group-specific ccREs
lop.analysis[[3]] <- my.function3(folder="analysis")
lop.analysis.all[[3]] <- my.function3(folder = "analysis.all")

# 4. number of ccREs intersecting each intron 
lop.analysis[[4]] <- my.function4(folder="analysis")
lop.analysis.all[[4]] <- my.function4(folder="analysis.all")

# 5. length of introns overlapping group-specific ccREs
lop.analysis[[5]] <- my.function5(folder="analysis")
lop.analysis.all[[5]] <- my.function5(folder="analysis.all")

pdf("~/public_html/enhancers_neural_development/plots/fig.12.analysis.pdf", width=8, height=4.5)
for ( i in 1:5) {
  
  print (lop.analysis[[i]])
}
dev.off()

pdf("~/public_html/enhancers_neural_development/plots/fig.12.analysis.all.pdf", width=8, height=4.5)
for ( i in 1:5) {
  
  print (lop.analysis.all[[i]])
}
dev.off()

