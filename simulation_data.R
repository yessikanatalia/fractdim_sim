### Libraries ###
library(dplyr)
library(tidyr)
library(fractaldim)
library(zoo)
library(forecast)
library(cowplot)
library(ggplot2)
library(scales)
library(ggforce)

### Set directory ###
dir <- ""

### Function for fractal dimension ###
fract.d <- function (data, window, method) {
  rollapply(data, width = window, 
            FUN = function (x) {
              return(tryCatch(fd.estimate(x,method=method, keep.loglog=TRUE,
                                          plot.loglog=FALSE, plot.allpoints=TRUE,
                                          nlags="auto")$fd, 
                              error=function(e) NA))},
            align="right", partial=T)
  
}

fillNA <- function(x) {
  ifelse(is.na(x),1,x)
}

### Generate initial incidence rate ##
lamb_pois <- seq(0.01,0.31,0.02) # lambda for the Poisson distribution
Day <- c(1:365) # length of time series data
Methods <- c("boxcount", "hallwood", "variogram", "madogram") # fractal dimension estimator
Window <- c(7,14,21) # sliding window
n_list <- list()
p_list <- list()
for (i in 1:length (lamb_pois))
{
  set.seed(0)
  n0 = arima.sim(list(order = c(0,0,0)), 
                 rand.gen = function(x) rpois(x, lambda = lamb_pois[i]),
                 n = length(Day))
  nbase = data.frame(cbind(Day,n0))
  n_list[[i]] = data.frame(crossing(Day, Methods, Window))
  n_list[[i]]$n0 = nbase$n0[match(n_list[[i]]$Day,nbase$Day)]
  n_list[[i]]$sim = i
  # Calculate fractal dimension of each time series data
  for (j in seq_along(Methods))
  {
    for (k in seq_along(Window))
    {
      subset = n_list[[i]]%>%
        filter(Methods==Methods[j],
               Window==Window[k])%>%
        arrange(Day)
      n_list[[i]]$fd[n_list[[i]]$Methods==Methods[j]&n_list[[i]]$Window==Window[k]] = 
        fillNA(fract.d(subset$n0, window = as.numeric(Window[k]),method = Methods[j]))
    }
  }
  
  # Text to annotate graphs
  graphLabels = n_list[[i]]%>%
    group_by(Methods,Window)%>%
    summarise(mean = round(mean(fd),2),
              variance = round(var(fd),2),
              acf = round(acf(fd,lag.max = 1, plot = FALSE,
                              na.action = na.pass)$acf[2],2) )
  
  # Create plot of incidence rate and loess regression of local fractal dimension
  p.7 <- ggplot(n_list[[i]][n_list[[i]]$Window==7,],aes(x=Day)) + 
    geom_smooth(aes(y=fd,color=Methods), method = "loess", 
                span=0.04, size=1, se = FALSE, method.args = list(degree=1)) +
    scale_color_manual(values = c("boxcount" = "red",
                                  "hallwood" = "forestgreen",
                                  "variogram" = "purple",
                                  "madogram" = "blue")) +
    theme_bw() +
    labs(title = paste("Sliding window = 7 days")) +
    scale_x_continuous(name="") +
    scale_y_continuous("Fractal dimension", 
                       limits=c(1,2)) +
    theme(axis.text.y=element_text(color=c("black")),
          axis.ticks.y=element_line(color=c("black")),
          axis.title.y=element_text(color=c("black")),
          axis.text.x = element_blank(),
          text = element_text(size = 10),
          plot.title = element_text(color = "black", size = 12, face = "bold"),
          legend.position = "none") 
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 2,label=paste0("Mean = ",mean)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 1.8,label=paste0("Var = ",variance)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 1.6,label=paste0("ACF = ",acf)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 2,label=paste0("Mean = ",mean)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 1.8,label=paste0("Var = ",variance)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 1.6,label=paste0("ACF = ",acf)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 2,label=paste0("Mean = ",mean)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 1.8,label=paste0("Var = ",variance)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 1.6,label=paste0("ACF = ",acf)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 2,label=paste0("Mean = ",mean)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 1.8,label=paste0("Var = ",variance)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 1.6,label=paste0("ACF = ",acf)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1)
  
  p.14 <- ggplot(n_list[[i]][n_list[[i]]$Window==14,],aes(x=Day)) + 
    geom_smooth(aes(y=fd,color=Methods), method = "loess", 
                span=0.04, size=1, se = FALSE, method.args = list(degree=1)) +
    scale_color_manual(values = c("boxcount" = "red",
                                  "hallwood" = "forestgreen",
                                  "variogram" = "purple",
                                  "madogram" = "blue")) +
    theme_bw() +
    labs(title = paste("Sliding window = 14 days")) +
    scale_x_continuous(name="") +
    scale_y_continuous("Fractal dimension", 
                       limits=c(1,2)) +
    theme(axis.text.y=element_text(color=c("black")),
          axis.ticks.y=element_line(color=c("black")),
          axis.title.y=element_text(color=c("black")),
          axis.text.x = element_blank(),
          text = element_text(size = 10),
          plot.title = element_text(color = "black", size = 12, face = "bold"),
          legend.position = "none") 
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 2,label=paste0("Mean = ",mean)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 1.8,label=paste0("Var = ",variance)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 1.6,label=paste0("ACF = ",acf)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 2,label=paste0("Mean = ",mean)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 1.8,label=paste0("Var = ",variance)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 1.6,label=paste0("ACF = ",acf)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 2,label=paste0("Mean = ",mean)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 1.8,label=paste0("Var = ",variance)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 1.6,label=paste0("ACF = ",acf)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 2,label=paste0("Mean = ",mean)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 1.8,label=paste0("Var = ",variance)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 1.6,label=paste0("ACF = ",acf)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1)
  
  p.21 <- ggplot(n_list[[i]][n_list[[i]]$Window==21,],aes(x=Day)) + 
    geom_smooth(aes(y=fd,color=Methods), method = "loess", 
                span=0.04, size=1, se = FALSE, method.args = list(degree=1)) +
    scale_color_manual(values = c("boxcount" = "red",
                                  "hallwood" = "forestgreen",
                                  "variogram" = "purple",
                                  "madogram" = "blue")) +
    theme_bw() +
    labs(title = paste("Sliding window = 21 days")) +
    scale_x_continuous(name="") +
    scale_y_continuous("Fractal dimension", 
                       limits=c(1,2)) +
    theme(axis.text.y=element_text(color=c("black")),
          axis.ticks.y=element_line(color=c("black")),
          axis.title.y=element_text(color=c("black")),
          axis.text.x = element_blank(),
          text = element_text(size = 10),
          plot.title = element_text(color = "black", size = 12, face = "bold"),
          legend.position = "none") 
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 2,label=paste0("Mean = ",mean)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 1.8,label=paste0("Var = ",variance)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="boxcount",], 
    #           aes(x = 100, y = 1.6,label=paste0("ACF = ",acf)), color="red" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 2,label=paste0("Mean = ",mean)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 1.8,label=paste0("Var = ",variance)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="hallwood",], 
    #           aes(x = 180, y = 1.6,label=paste0("ACF = ",acf)), color="forestgreen" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 2,label=paste0("Mean = ",mean)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 1.8,label=paste0("Var = ",variance)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="madogram",], 
    #           aes(x = 260, y = 1.6,label=paste0("ACF = ",acf)), color="blue" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 2,label=paste0("Mean = ",mean)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 1.8,label=paste0("Var = ",variance)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1) +
    # geom_text(data = graphLabels[graphLabels$Methods=="variogram",], 
    #           aes(x = 340, y = 1.6,label=paste0("ACF = ",acf)), color="purple" ,
    #           lineheight = 0.75, size = 3, hjust=1)
  
  p.cases <- ggplot(n_list[[i]][n_list[[i]]$Window==7,],aes(x=Day)) + 
    geom_line(aes(y=n0), size=1, color = "black")+
    theme_bw() +
    scale_x_continuous(name="Days") +
    scale_y_continuous("Incidence rate", 
                       limits=c(0,3)) +
    theme(axis.text.y=element_text(color=c("black")),
          axis.ticks.y=element_line(color=c("black")),
          axis.title.y=element_text(color=c("black")),
          axis.text.x = element_text(angle = 0, hjust = 0, color = "black"),
          text = element_text(size = 10),
          plot.title = element_text(color = "black", size = 12, face = "bold"),
          legend.position = "none")
  
  title <- ggdraw() + draw_label(paste("\u03BB","=",lamb_pois[i]), fontface='bold')
  
  p_list[[i]] <- cowplot::plot_grid(title,p.7,p.14,p.21,p.cases,
                                    ncol=1, align = "v",
                                    rel_heights=c(0.2,1,1,1,1))
}

ggsave(cowplot::plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],
                          p_list[[4]],p_list[[5]],p_list[[6]],
                          p_list[[7]],p_list[[8]],p_list[[9]],
                          p_list[[10]],p_list[[11]],p_list[[12]],
                          p_list[[13]],p_list[[14]],p_list[[15]],
                          p_list[[16]],
                          ncol = 4),
       file = paste0(dir, "individual_curves.png"), 
       width = 8000, height = 8000, units="px", device='png', dpi=300)

### Replication function ###
cases_sim <- function (N, Lambda, Methods, Window)
{
  nbase <- list()
  for (i in seq_along(Lambda))
  {
    n0 = arima.sim(list(order = c(0,0,0)), 
                   rand.gen = function(x) rpois(x, lambda = Lambda[i]),
                   n = N)
    data0 = data.frame(cbind(Day=c(1:N),n0))
    nbase[[i]] = data.frame(crossing(Day=c(1:N), Methods, Window, Lambda[i]))
    nbase[[i]]$n0 = data0$n0[match(nbase[[i]]$Day,data0$Day)]
    for (j in 1:length(Methods))
    {
      for (k in 1:length(Window))
      {
        subset = nbase[[i]]%>%
          filter(Methods==Methods[j],
                 Window==Window[k])%>%
          arrange(Day)
        nbase[[i]]$fd[nbase[[i]]$Methods==Methods[j]&nbase[[i]]$Window==Window[k]] = 
          fillNA(fract.d(subset$n0, window = as.numeric(Window[k]),
                         method = Methods[j]))
      }
    }
  }
  nbase_long <- bind_rows(nbase)
  names(nbase_long)[names(nbase_long) == "Lambda.i."] <- "Lambda"
  
  sumbase <- nbase_long%>%
    group_by(Methods,Window,Lambda)%>%
    summarise(mean.fd = round(mean(fd, na.rm=TRUE),2),
              variance.fd = round(var(fd),2),
              acf.fd = round(acf(fd,lag.max = 1, plot = FALSE,
                              na.action = na.pass)$acf[2],2))
  sumbase
}

### Create initial dataset with 1000 replication ###
N <- 365
lamb_pois <- seq(0.01,0.31,0.02)
methods <- c("boxcount", "hallwood", "variogram", "madogram")
window <- c(7,14,21)

sim_list <- list()
for (i in 1:1000)
{
   sim_list[[i]] <- cases_sim(N=N,Lambda = lamb_pois, Methods = methods, Window = window)
   sim_list[[i]]$sim <- i
}

sim_long <- bind_rows(sim_list)%>%
  mutate(meth.class = case_when(Methods== "boxcount" ~ "Boxcount",
                                Methods== "hallwood" ~ "Hall-Wood",
                                Methods== "variogram" ~ "Variogram",
                                Methods== "madogram" ~ "Madogram"))
                     
### Save initial dataset for further analysis ### 
write.csv(sim_long,paste0(dir,"sim_long.csv"), row.names = FALSE)
