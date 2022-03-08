### Libraries
library(dplyr)
library(tidyr)
library(fractaldim)
library(zoo)
library(forecast)
library(cowplot)
library(ggplot2)
library(scales)
library(ggforce)
library(factoextra)

### Function for fractal dimension
fract.d = function (data, window, method) {
  rollapply(data, width = window, 
            FUN = function (x) {
              return(tryCatch(fd.estimate(x,method=method, keep.loglog=TRUE,
                                          plot.loglog=FALSE, plot.allpoints=TRUE,
                                          nlags="auto")$fd, 
                              error=function(e) NA))},
            align="right", partial=T)
  
}

mse.d = function (data, window, method) {
  rollapply(data, width = window, 
            FUN = function (x) {
              sqrt(return(tryCatch(fd.estimate(x,method=method, keep.loglog=TRUE,
                                          plot.loglog=FALSE, plot.allpoints=TRUE,
                                          nlags="auto")$loglog[[1]][[1]][["lsq"]], 
                              error=function(e) NA)))},
            align="right", partial=T)
  
}

fillNA = function(x) {
  ifelse(is.na(x),1,x)
}

qfun <- function(x, q_1=0.25, q_2=0.5, q_3=0.75){
  q <- sort(c(quantile(x, q_1), quantile(x, q_2), quantile(x, q_3)))
  return(c(min(x), q[1], q[2], q[3], max(x)))
}

### Example of Poisson + MA model##
lamb_pois <- seq(0.01,0.5,0.02) #lambda for the Poisson distribution
n <- 365 #length of series
Day=c(1:365)
Methods=c("boxcount", "hallwood", "variogram", "madogram")
Window=c(7,14,21)

n_list=list()
p_list=list()
for (i in 1:18)
{
  set.seed(0)
  n0 = arima.sim(list(order = c(0,0,0)), 
                 rand.gen = function(x) rpois(x, lambda = lamb_pois[i]),
                 n = n)
  nbase = data.frame(cbind(Day,n0))
  n_list[[i]] = data.frame(crossing(Day, Methods, Window))
  n_list[[i]]$n0 = nbase$n0[match(n_list[[i]]$Day,nbase$Day)]
  n_list[[i]]$sim = i
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
      # n_list[[i]]$mse[n_list[[i]]$Methods==Methods[j]&n_list[[i]]$Window==Window[k]] = 
      #   mse.d(subset$n0, window = as.numeric(Window[k]),method = Methods[j])
    }
  }
  
  # text to annotate graphs
  graphLabels = n_list[[i]]%>%
    group_by(Methods,Window)%>%
    summarise(mean = round(mean(fd),2),
              variance = round(var(fd),2),
              acf = round(acf(fd,lag.max = 1, plot = FALSE,
                              na.action = na.pass)$acf[2],2) )
  
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
  
  #p.mse = ggplot(n_list[[i]],aes(x=day)) + 
  #  geom_point(aes(y=mse,colour=methods),size=1) +
  #  theme_bw() +
  #  labs(title = paste0("Lambda = ",lambda[i])) +
  #  scale_x_continuous(name="Day") +
  #  scale_y_continuous("RMSE",limits = c(0,10^-30)) +
  #  theme(axis.text.y.right=element_text(color=c("black")),
  #        axis.ticks.y.right=element_line(color=c("black")),
  #        axis.title.y.right=element_text(color=c("black")),
  #        axis.text.y=element_text(color=c("black")),
  #        axis.ticks.y=element_line(color=c("black")),
  #        axis.title.y=element_text(color=c("black")),
  #        axis.text.x = element_text(angle = 0, hjust = 0, color = "black"),
  #        text = element_text(size = 10),
  #        plot.title = element_text(color = "black", size = 12, face = "bold"),
  #        legend.key.height= unit(0.5, 'cm'),
  #        legend.key.width= unit(0.5, 'cm'),
  #        legend.key.size = unit(3, 'cm'),
  #        legend.title = element_text(size = 10),
  #        legend.text = element_text(size = 10)) +
  #  guides(color = guide_legend(override.aes = list(size = 1.5)))+
  #  facet_wrap(~Window, nrow=3)
}

ggsave(cowplot::plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],
                          p_list[[4]],p_list[[5]],p_list[[6]],
                          p_list[[7]],p_list[[8]],p_list[[9]],
                          p_list[[10]],p_list[[11]],p_list[[12]],
                          p_list[[13]],p_list[[14]],p_list[[15]],
                          p_list[[16]],
                          ncol = 4),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/indiv.png"), 
       width = 8000, height = 8000, units="px", device='png', dpi=300)

### Replication function
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
        # nbase$mse[nbase$Methods==Methods[j]&nbase$Window==Window[k]]] = 
        #   mse.d(subset$n0, window = as.numeric(Window[k]),method = Methods[j])
        
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
              # mean.mse = mean(mse, na.rm=TRUE))
  sumbase
}

### Create original dataset with 1000 replication
N=365
lamb_pois=seq(0.01,0.31,0.02)
methods=c("boxcount", "hallwood", "variogram", "madogram")
window=c(7,14,21)

# sim_list <- list()
# for (i in 1:1000)
# {
#   sim_list[[i]] <- cases_sim(N=N,Lambda = lamb_pois, Methods = methods, Window = window)
#   sim_list[[i]]$sim <- i
# }
# 
# sim_long <- bind_rows(sim_list)%>%
#   mutate(meth.class = case_when(Methods== "boxcount" ~ "Boxcount",
#                                 Methods== "hallwood" ~ "Hall-Wood",
#                                 Methods== "variogram" ~ "Variogram",
#                                 Methods== "madogram" ~ "Madogram"))
# 
# write.csv(sim_long,"G:/My Drive/COVID-19-Fractal dimension/data/sim_long_2.csv", row.names = FALSE)

sim_long <- read.csv("G:/My Drive/COVID-19-Fractal dimension/data/sim_long_2.csv")

p_list <- list()
for (i in 1:length(lamb_pois))
{
  data0 <- sim_long%>%
    filter (Lambda == 0.29)
  
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
  
  #p.mse.fd = ggplot(sim_long) +
  #  geom_boxplot(aes(x=Methods,y=mean.mse,fill=Methods)) +
  #  theme_bw() +
  #  labs(title = paste0("Lambda = ",Lambda[i])) +
  #  scale_y_continuous("Mean RMSE",limits = c(0,5*10^-31),
  #                     labels = scientific_format()) +
  #  theme(legend.position = "none") +
  #  facet_wrap(~Window, nrow=3, scales = "free_x")
  
  title <- ggdraw() + draw_label(paste("\u03BB","=",lamb_pois[15]), 
                                 fontface='bold')
  
  p_list[[15]] <- cowplot::plot_grid(title, p.mean,p.var,p.acf,
                                    ncol=1, align = "v",
                                    rel_heights=c(0.15,1,1,1))
}

ggsave(cowplot::plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],
                          ncol=4),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/rep.1.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

ggsave(cowplot::plot_grid(p_list[[5]],p_list[[6]],p_list[[7]],p_list[[8]],
                          ncol=4),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/rep.2.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

ggsave(cowplot::plot_grid(p_list[[9]],p_list[[10]],p_list[[11]],p_list[[12]],
                          ncol=4),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/rep.3.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)

ggsave(cowplot::plot_grid(p_list[[13]],p_list[[14]],p_list[[15]],p_list[[16]],
                          ncol=4),
       file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/rep.4.png"), 
       width = 8000, height = 6000, units="px", device='png', dpi=300)



### Cluster analysis (k-means)
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
    
    # ggsave(plot_grid(fig.7.all,fig.14.all,fig.21.all,ncol=2),
    #        file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/sim.kmean.",i,
    #                      methods[j],".png"),
    #        width = 4000, height = 3000, units="px", device='png', dpi=300)
    
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
           file = paste0("G:/My Drive/COVID-19-Fractal dimension/figures/sim.kmean.rev.",i,
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

sim_long_agg <- crossing(sim=c(1:1000), methods, window, cluster=c(1:4),
                         mean.centroid=NA, var.centroid=NA, acf.centroid=NA,
                         mean.mean=NA, var.mean=NA, acf.mean=NA,
                         mean.stat=NA, var.stat=NA, acf.stat=NA)
for (i in 1:1000)
{
  for (j in seq_along(methods))
  {
    for(k in seq_along(window))
    {
      data.0 <- sim_long%>%
        filter(sim==i,
               Methods==methods[j],
               Window==window[k])
      
      data.0 <- data.frame(na.omit(data.0[,c(3:6)]))
      rownames(data.0) <- data.0[,1]
      data.0[,1] <- NULL
      
      # Calculate the clusters
      clust.0 <- kmeans(data.0, centers = 4, nstart = 25)
      clust.0.center <- data.frame(clust.0$centers)
      colnames(clust.0.center) <- c("mean.centroid","var.centroid","acf.centroid")
      clust.0.center <- tibble::rownames_to_column(clust.0.center, "cluster")
      clust.0.center$cluster <- as.numeric(clust.0.center$cluster)
      
      # Fill the data
      sim_long_agg$mean.centroid[sim_long_agg$sim==i&sim_long_agg$methods==methods[j]&
                               sim_long_agg$window==window[k]] <- clust.0.center$mean.centroid[match(sim_long_agg$cluster,clust.0.center$cluster)]
      sim_long_agg$var.centroid[sim_long_agg$sim==i&sim_long_agg$methods==methods[j]&
                               sim_long_agg$window==window[k]] <- clust.0.center$var.centroid[match(sim_long_agg$cluster,clust.0.center$cluster)]
      sim_long_agg$acf.centroid[sim_long_agg$sim==i&sim_long_agg$methods==methods[j]&
                               sim_long_agg$window==window[k]] <- clust.0.center$acf.centroid[match(sim_long_agg$cluster,clust.0.center$cluster)]
      sim_long_agg$mean.mean[sim_long_agg$sim==i&sim_long_agg$methods==methods[j]&
                                  sim_long_agg$window==window[k]] <- mean(as.numeric(data.0[,1]))
      sim_long_agg$var.mean[sim_long_agg$sim==i&sim_long_agg$methods==methods[j]&
                               sim_long_agg$window==window[k]] <- mean(as.numeric(data.0[,2]))
      sim_long_agg$acf.mean[sim_long_agg$sim==i&sim_long_agg$methods==methods[j]&
                               sim_long_agg$window==window[k]] <- mean(as.numeric(data.0[,3]))
    }
  }
}

sim_long_agg$mean.stat[sim_long_agg$mean.centroid >= sim_long_agg$mean.mean] <- 1
sim_long_agg$mean.stat[sim_long_agg$mean.centroid < sim_long_agg$mean.mean] <- 0
sim_long_agg$var.stat[sim_long_agg$var.centroid >= sim_long_agg$var.mean] <- 1
sim_long_agg$var.stat[sim_long_agg$var.centroid < sim_long_agg$var.mean] <- 0
sim_long_agg$acf.stat[sim_long_agg$acf.centroid >= sim_long_agg$acf.mean] <- 1
sim_long_agg$acf.stat[sim_long_agg$acf.centroid < sim_long_agg$acf.mean] <- 0
