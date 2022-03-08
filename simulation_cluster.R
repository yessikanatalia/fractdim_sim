### Libraries ###
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(scales)
library(ggforce)
library(factoextra)

### Set directory ###
dir <- ""

### Load replication data ###
sim_long <- read.csv(paste0(dir,"sim_long.csv"))
lamb_pois <- seq(0.01,0.31,0.02)
methods <- c("boxcount", "hallwood", "variogram", "madogram")
window <- c(7,14,21)

### Generate boxplots ###
p_list <- list()
for (i in 1:length(lamb_pois))
{
  data0 <- sim_long%>%
    filter (Lambda == lamb_pois[i])
  
  p.mean = ggplot(data0) +
    geom_boxplot(aes(x=meth.class,y=mean.fd,fill=meth.class)) +
    scale_fill_manual(values = c("Boxcount" = "red",
                                  "Hall-Wood" = "forestgreen",
                                  "Variogram" = "purple",
                                  "Madogram" = "blue")) +
    theme_bw() +
    scale_y_continuous("Mean FD", limits = c(1,2)) +
    labs(x = "") +
    theme(legend.position = "none") +
    facet_wrap(~Window, nrow=3)
  
  p.var = ggplot(data0) +
    geom_boxplot(aes(x=meth.class,y=variance.fd,fill=meth.class)) +
    scale_fill_manual(values = c("Boxcount" = "red",
                                 "Hall-Wood" = "forestgreen",
                                 "Variogram" = "purple",
                                 "Madogram" = "blue")) +
    theme_bw() +
    scale_y_continuous("Variance FD", limits = c(0,0.3)) +
    labs(x = "") +
    theme(legend.position = "none") +
    facet_wrap(~Window, nrow=3)
  
  p.acf = ggplot(data0) +
    geom_boxplot(aes(x=meth.class,y=acf.fd,fill=meth.class)) +
    scale_fill_manual(values = c("Boxcount" = "red",
                                 "Hall-Wood" = "forestgreen",
                                 "Variogram" = "purple",
                                 "Madogram" = "blue")) +
    theme_bw() +
    scale_y_continuous("ACF FD", limits = c(-1,1)) +
    labs(x = "") +
    theme(legend.position = "none") +
    facet_wrap(~Window, nrow=3)
  
  title <- ggdraw() + draw_label(paste("\u03BB","=",lamb_pois[15]), 
                                 fontface='bold')
  
  p_list[[i]] <- cowplot::plot_grid(title, p.mean,p.var,p.acf,
                                    ncol=1, align = "v",
                                    rel_heights=c(0.15,1,1,1))
}

ggsave(cowplot::plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],
                          ncol=4),
       file = paste0(dir,"rep.1.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

ggsave(cowplot::plot_grid(p_list[[5]],p_list[[6]],p_list[[7]],p_list[[8]],
                          ncol=4),
       file = paste0(dir,"rep.2.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

ggsave(cowplot::plot_grid(p_list[[9]],p_list[[10]],p_list[[11]],p_list[[12]],
                          ncol=4),
       file = paste0(dir,"rep.3.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

ggsave(cowplot::plot_grid(p_list[[13]],p_list[[14]],p_list[[15]],p_list[[16]],
                          ncol=4),
       file = paste0(dir,"rep.4.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

### Cluster analysis (k-means) ###
## Check for optimal number of cluster + save figures
for (i in 1:5) # choose one or several simulations to save time
{
  for (j in seq_along(methods))
  {
    subset <- sim_long%>%
      filter(Methods==methods[j])
    
    # Subset data per sliding window
    data.7 <- subset%>%
      filter(Window==7)
    data.7 <- data.frame(na.omit(data.7[,c(4:6)]))
    
    data.14 <- subset%>%
      filter(Window==14)
    data.14 <- data.frame(na.omit(data.14[,c(4:6)]))
    
    data.21 <- subset%>%
      filter(Window==21)
    data.21 <- data.frame(na.omit(data.21[,c(4:6)]))
    
    # Optimal number of cluster (elbow method)
    print(fviz_nbclust(data.7, kmeans, method = "wss") +
            labs(subtitle = paste0(subset$meth.class, ", sliding window = 7 days")))
    print(fviz_nbclust(data.14, kmeans, method = "wss") +
            labs(subtitle = paste0(subset$meth.class, ", sliding window = 14 days")))
    print(fviz_nbclust(data.21, kmeans, method = "wss") +
            labs(subtitle = paste0(subset$meth.class, ", sliding window = 21 days")))
  }
}

## Create cluster figures
for (i in 1:5) # choose one or several simulations to save time
{
  for (j in seq_along(methods))
  {
    subset <- sim_long%>%
      filter(sim==i,
             Methods==methods[j])
    
    # Subset data per sliding window
    data.7 <- subset%>%
      filter(Window==7)
    data.7 <- data.frame(na.omit(data.7[,c(3:6)]))
    rownames(data.7) <- data.7[,1]
    data.7[,1] <- NULL
    
    data.14 <- subset%>%
      filter(Window==14)
    data.14 <- data.frame(na.omit(data.14[,c(3:6)]))
    rownames(data.14) <- data.14[,1]
    data.14[,1] <- NULL
    
    data.21 <- subset%>%
      filter(Window==21)
    data.21 <- data.frame(na.omit(data.21[,c(3:6)]))
    rownames(data.21) <- data.21[,1]
    data.21[,1] <- NULL
    
    # Calculate the clusters
    clust.7 <- kmeans(data.7, centers = 4, nstart = 25)
    clust.14 <- kmeans(data.14, centers = 4, nstart = 25)
    clust.21 <- kmeans(data.21, centers = 4, nstart = 25)
    
    # Create initial figures
    fig.7 <- fviz_cluster(clust.7, data = data.7, geom = c("point")) +
      scale_color_brewer('Cluster', palette='Set1') + 
      scale_fill_brewer('Cluster', palette='Set1') +
      scale_shape_manual('Cluster', values=c(22,23,24,25)) + 
      labs(title = paste0(subset$meth.class, ", sliding window = 7 days")) +
      theme_bw()
    fig.7.all <- fig.7 + 
      geom_text(data = fig.7$data, aes(x=x, y=y, label=name, colour=cluster),
                size = 3, vjust=-1, show.legend = F)
    
    fig.14 <- fviz_cluster(clust.14, data = data.14, geom = c("point")) +
      scale_color_brewer('Cluster', palette='Set1') + 
      scale_fill_brewer('Cluster', palette='Set1') +
      scale_shape_manual('Cluster', values=c(22,23,24,25)) + 
      labs(title = paste0(subset$meth.class, ", sliding window = 14 days")) +
      theme_bw()
    fig.14.all <- fig.14 + 
      geom_text(data = fig.14$data, aes(x=x, y=y, label=name, colour=cluster),
                size = 3, vjust=-1, show.legend = F)
    
    fig.21 <- fviz_cluster(clust.21, data = data.21, geom = c("point")) +
      scale_color_brewer('Cluster', palette='Set1') + 
      scale_fill_brewer('Cluster', palette='Set1') +
      scale_shape_manual('Cluster', values=c(22,23,24,25)) + 
      labs(title = paste0(subset$meth.class, ", sliding window = 21 days")) +
      theme_bw()
    fig.21.all <- fig.21 + 
      geom_text(data = fig.21$data, aes(x=x, y=y, label=name, colour=cluster),
                size = 3, vjust=-1, show.legend = F)
    
    # Extract the clusters & centroids
    clust.7.center <- data.frame(clust.7$centers)
    colnames(clust.7.center) <- c("mean.centroid","var.centroid","acf.centroid")
    clust.7.center <- tibble::rownames_to_column(clust.7.center, "cluster")
    clust.7.center$cluster <- as.numeric(clust.7.center$cluster)
    
    clusters.7 <- data.frame(clust.7$cluster)
    colnames(clusters.7) <- c("cluster")
    clusters.7 <- tibble::rownames_to_column(clusters.7, "Lambda")
    
    clust.14.center <- data.frame(clust.14$centers)
    colnames(clust.14.center) <- c("mean.centroid","var.centroid","acf.centroid")
    clust.14.center <- tibble::rownames_to_column(clust.14.center, "cluster")
    clust.14.center$cluster <- as.numeric(clust.14.center$cluster)
    
    clusters.14 <- data.frame(clust.14$cluster)
    colnames(clusters.14) <- c("cluster")
    clusters.14 <- tibble::rownames_to_column(clusters.14, "Lambda")
    
    clust.21.center <- data.frame(clust.21$centers)
    colnames(clust.21.center) <- c("mean.centroid","var.centroid","acf.centroid")
    clust.21.center <- tibble::rownames_to_column(clust.21.center, "cluster")
    clust.21.center$cluster <- as.numeric(clust.21.center$cluster)
    
    clusters.21 <- data.frame(clust.21$cluster)
    colnames(clusters.21) <- c("cluster")
    clusters.21 <- tibble::rownames_to_column(clusters.21, "Lambda")
    
    # Combine with original data
    clusters.7.all <- left_join(clusters.7, clust.7.center, by=c("cluster"))
    clusters.7.all$mean.mean <- mean(as.numeric(data.7[,1]))
    clusters.7.all$mean.var <- mean(as.numeric(data.7[,2]))
    clusters.7.all$mean.acf <- mean(as.numeric(data.7[,3]))
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "Low mean, low variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "Low mean, low variance, high acf"
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "Low mean, high variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "Low mean, high variance, high acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "High mean, low variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "High mean, low variance, high acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "High mean, high variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "High mean, high variance, high acf"
    
    clusters.14.all <- left_join(clusters.14, clust.14.center, by=c("cluster"))
    clusters.14.all$mean.mean <- mean(as.numeric(data.14[,1]))
    clusters.14.all$mean.var <- mean(as.numeric(data.14[,2]))
    clusters.14.all$mean.acf <- mean(as.numeric(data.14[,3]))
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid < clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "Low mean, low variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid < clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "Low mean, low variance, high acf"
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "Low mean, high variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "Low mean, high variance, high acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid < clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "High mean, low variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid < clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "High mean, low variance, high acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "High mean, high variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                            clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                            clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "High mean, high variance, high acf"
    
    clusters.21.all <- left_join(clusters.21, clust.21.center, by=c("cluster"))
    clusters.21.all$mean.mean <- mean(as.numeric(data.21[,1]))
    clusters.21.all$mean.var <- mean(as.numeric(data.21[,2]))
    clusters.21.all$mean.acf <- mean(as.numeric(data.21[,3]))
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid < clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "Low mean, low variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid < clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "Low mean, low variance, high acf"
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "Low mean, high variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "Low mean, high variance, high acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid < clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "High mean, low variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid < clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "High mean, low variance, high acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "High mean, high variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                            clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                            clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "High mean, high variance, high acf"
    
    # Create revised figures
    fig.7.label=c(unique(clusters.7.all$label[clusters.7.all$cluster==1]),
                  unique(clusters.7.all$label[clusters.7.all$cluster==2]),
                  unique(clusters.7.all$label[clusters.7.all$cluster==3]),
                  unique(clusters.7.all$label[clusters.7.all$cluster==4]))
    
    fig.7.edit <- fviz_cluster(clust.7, data = data.7, geom = c("point")) +
      scale_color_brewer('Cluster', palette='Set1',
                         labels=fig.7.label) + 
      scale_fill_brewer('Cluster', palette='Set1',
                        labels=fig.7.label) +
      scale_shape_manual('Cluster', values=c(22,23,24,25),
                         labels=fig.7.label) +
      labs(title = paste0(subset$meth.class, ", sliding window = 7 days")) +
      theme_bw()
    
    fig.7.rev <- fig.7.edit +
      geom_text(data = fig.7$data, aes(x=x, y=y, label=name, colour=cluster),
                size = 3, vjust=-1, show.legend = F)
    
    fig.14.label=c(unique(clusters.14.all$label[clusters.14.all$cluster==1]),
                   unique(clusters.14.all$label[clusters.14.all$cluster==2]),
                   unique(clusters.14.all$label[clusters.14.all$cluster==3]),
                   unique(clusters.14.all$label[clusters.14.all$cluster==4]))
    
    fig.14.edit <- fviz_cluster(clust.14, data = data.14, geom = c("point")) +
      scale_color_brewer('Cluster', palette='Set1',
                         labels=fig.14.label) + 
      scale_fill_brewer('Cluster', palette='Set1',
                        labels=fig.14.label) +
      scale_shape_manual('Cluster', values=c(22,23,24,25),
                         labels=fig.14.label) +
      labs(title = paste0(subset$meth.class, ", sliding window = 14 days")) +
      theme_bw() 
    
    fig.14.rev <- fig.14.edit + 
      geom_text(data = fig.14$data, aes(x=x, y=y, label=name, colour=cluster),
                size = 3, vjust=-1, show.legend = F)
    
    fig.21.label=c(unique(clusters.21.all$label[clusters.21.all$cluster==1]),
                   unique(clusters.21.all$label[clusters.21.all$cluster==2]),
                   unique(clusters.21.all$label[clusters.21.all$cluster==3]),
                   unique(clusters.21.all$label[clusters.21.all$cluster==4]))
    
    fig.21.edit <- fviz_cluster(clust.21, data = data.21, geom = c("point")) +
      scale_color_brewer('Cluster', palette='Set1',
                         labels=fig.21.label) + 
      scale_fill_brewer('Cluster', palette='Set1',
                        labels=fig.21.label) +
      scale_shape_manual('Cluster', values=c(22,23,24,25),
                         labels=fig.21.label) +
      labs(title = paste0(subset$meth.class, ", sliding window = 21 days")) +
      theme_bw() 
    fig.21.rev <- fig.21.edit +
      geom_text(data = fig.21$data, aes(x=x, y=y, label=name, colour=cluster),
                size = 3, vjust=-1, show.legend = F)
    
    ggsave(plot_grid(fig.7.rev,fig.14.rev,fig.21.rev,ncol=2),
           file = paste0(dir,"sim.kmean.rev.",i,".",
                         methods[j],".png"),
           width = 4000, height = 3000, units="px", device='png', dpi=300)
    
  }
}

## Validity check
valid_check <- crossing(sim=c(1:1000),methods,window,check=NA)

for (i in 1:1000)
{
  for (j in seq_along(methods))
  {
    subset <- sim_long%>%
      filter(sim==i,
             Methods==methods[j])
    
    # Subset data per sliding window
    data.7 <- subset%>%
      filter(Window==7)
    data.7 <- data.frame(na.omit(data.7[,c(3:6)]))
    rownames(data.7) <- data.7[,1]
    data.7[,1] <- NULL
    
    data.14 <- subset%>%
      filter(Window==14)
    data.14 <- data.frame(na.omit(data.14[,c(3:6)]))
    rownames(data.14) <- data.14[,1]
    data.14[,1] <- NULL
    
    data.21 <- subset%>%
      filter(Window==21)
    data.21 <- data.frame(na.omit(data.21[,c(3:6)]))
    rownames(data.21) <- data.21[,1]
    data.21[,1] <- NULL
    
    # Calculate the clusters
    clust.7 <- kmeans(data.7, centers = 4, nstart = 25)
    clust.14 <- kmeans(data.14, centers = 4, nstart = 25)
    clust.21 <- kmeans(data.21, centers = 4, nstart = 25)
    
    # Extract the clusters & centroids
    clust.7.center <- data.frame(clust.7$centers)
    colnames(clust.7.center) <- c("mean.centroid","var.centroid","acf.centroid")
    clust.7.center <- tibble::rownames_to_column(clust.7.center, "cluster")
    clust.7.center$cluster <- as.numeric(clust.7.center$cluster)
    
    clusters.7 <- data.frame(clust.7$cluster)
    colnames(clusters.7) <- c("cluster")
    clusters.7 <- tibble::rownames_to_column(clusters.7, "Lambda")
    
    clust.14.center <- data.frame(clust.14$centers)
    colnames(clust.14.center) <- c("mean.centroid","var.centroid","acf.centroid")
    clust.14.center <- tibble::rownames_to_column(clust.14.center, "cluster")
    clust.14.center$cluster <- as.numeric(clust.14.center$cluster)
    
    clusters.14 <- data.frame(clust.14$cluster)
    colnames(clusters.14) <- c("cluster")
    clusters.14 <- tibble::rownames_to_column(clusters.14, "Lambda")
    
    clust.21.center <- data.frame(clust.21$centers)
    colnames(clust.21.center) <- c("mean.centroid","var.centroid","acf.centroid")
    clust.21.center <- tibble::rownames_to_column(clust.21.center, "cluster")
    clust.21.center$cluster <- as.numeric(clust.21.center$cluster)
    
    clusters.21 <- data.frame(clust.21$cluster)
    colnames(clusters.21) <- c("cluster")
    clusters.21 <- tibble::rownames_to_column(clusters.21, "Lambda")
    
    # Combine with original data
    clusters.7.all <- left_join(clusters.7, clust.7.center, by=c("cluster"))
    clusters.7.all$mean.mean <- mean(as.numeric(data.7[,1]))
    clusters.7.all$mean.var <- mean(as.numeric(data.7[,2]))
    clusters.7.all$mean.acf <- mean(as.numeric(data.7[,3]))
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "Low mean, low variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "Low mean, low variance, high acf"
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "Low mean, high variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid < clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "Low mean, high variance, high acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "High mean, low variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid < clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "High mean, low variance, high acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid < clusters.7.all$mean.acf] <- "High mean, high variance, low acf"
    clusters.7.all$label[clusters.7.all$mean.centroid >= clusters.7.all$mean.mean &
                           clusters.7.all$var.centroid >= clusters.7.all$mean.var &
                           clusters.7.all$acf.centroid >= clusters.7.all$mean.acf] <- "High mean, high variance, high acf"

    clusters.14.all <- left_join(clusters.14, clust.14.center, by=c("cluster"))
    clusters.14.all$mean.mean <- mean(as.numeric(data.14[,1]))
    clusters.14.all$mean.var <- mean(as.numeric(data.14[,2]))
    clusters.14.all$mean.acf <- mean(as.numeric(data.14[,3]))
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid < clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "Low mean, low variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid < clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "Low mean, low variance, high acf"
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "Low mean, high variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid < clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "Low mean, high variance, high acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid < clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "High mean, low variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid < clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "High mean, low variance, high acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid < clusters.14.all$mean.acf] <- "High mean, high variance, low acf"
    clusters.14.all$label[clusters.14.all$mean.centroid >= clusters.14.all$mean.mean &
                           clusters.14.all$var.centroid >= clusters.14.all$mean.var &
                           clusters.14.all$acf.centroid >= clusters.14.all$mean.acf] <- "High mean, high variance, high acf"
    
    clusters.21.all <- left_join(clusters.21, clust.21.center, by=c("cluster"))
    clusters.21.all$mean.mean <- mean(as.numeric(data.21[,1]))
    clusters.21.all$mean.var <- mean(as.numeric(data.21[,2]))
    clusters.21.all$mean.acf <- mean(as.numeric(data.21[,3]))
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid < clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "Low mean, low variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid < clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "Low mean, low variance, high acf"
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "Low mean, high variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid < clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "Low mean, high variance, high acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid < clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "High mean, low variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid < clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "High mean, low variance, high acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid < clusters.21.all$mean.acf] <- "High mean, high variance, low acf"
    clusters.21.all$label[clusters.21.all$mean.centroid >= clusters.21.all$mean.mean &
                           clusters.21.all$var.centroid >= clusters.21.all$mean.var &
                           clusters.21.all$acf.centroid >= clusters.21.all$mean.acf] <- "High mean, high variance, high acf"
    
    if (length(unique(clusters.7.all$cluster))==length(unique(clusters.7.all$label)))
    {valid_check$check[valid_check$sim==i&valid_check$methods==methods[j]&valid_check$window==7] <- 1} else
      {valid_check$check[valid_check$sim==i&valid_check$methods==methods[j]&valid_check$window==7] <-0}

    if (length(unique(clusters.14.all$cluster))==length(unique(clusters.14.all$label)))
    {valid_check$check[valid_check$sim==i&valid_check$methods==methods[j]&valid_check$window==14] <- 1} else
      {valid_check$check[valid_check$sim==i&valid_check$methods==methods[j]&valid_check$window==14] <-0}

    if (length(unique(clusters.21.all$cluster))==length(unique(clusters.21.all$label)))
    {valid_check$check[valid_check$sim==i&valid_check$methods==methods[j]&valid_check$window==21] <- 1} else
      {valid_check$check[valid_check$sim==i&valid_check$methods==methods[j]&valid_check$window==21] <-0}

  }  
}  

valid_check_agg <- valid_check%>%
  group_by(methods,window)%>%
  summarise(perc=length(check[check==1])/length(check)*100)
