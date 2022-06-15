Treatment_name =  c("Perch- / Nutrients+","Perch+ / Nutrients+","Perch- / Nutrients-","Perch+ / Nutrients-")

# ------------------------------
# EXTRACTION OF DATA FROM MODEL
# ------------------------------

# Extract phi from model (survival)
get_phi <- function(model){
  model <- as.mcmc(model)
  phi = data.frame()
  for (n in 1:2){
    for (tr in 1:4){
      for (gs in 1:3){
        for (e in 1:2){
          phi <- rbind(phi, 
                       data.frame(model[[n]][,paste("phi",as.character(gs),"[",as.character(tr),",",as.character(e),"]", sep = "")],
                                  Size = gs,
                                  Treatment = tr,
                                  Exp_time = e))
        }
      }
    }
  }
  
  return (phi %>% 
            group_by (Size,Treatment,Exp_time) %>% 
            summarize(mean = mean(var1), min = quantile(var1,0.025), max = quantile(var1,0.975)) %>%
            mutate(Size = factor(Size, levels = c("1","2","3"), labels = c("Small","Medium","Large")),
                   Treatment  = factor(Treatment, levels = c("1","2","3","4"), labels = Treatment_name),
                   Exp_time = factor(Exp_time, levels = c("1","2"), 
                                     labels = c("Before treatment","After treatment"))))
}

# Extract p from model (capture probability)
get_p <- function(model){
  model <- as.mcmc(model)
  p = data.frame()
  for (n in 1:2){
    for (gs in 1:3){
      p <- rbind(p, 
                 data.frame(model[[n]][,paste("p",as.character(gs), sep = "")],
                            Size = gs))
    }
  }
  
  return(p %>% 
           group_by (Size) %>% 
           summarize(mean = mean(var1), min = quantile(var1,0.025), max = quantile(var1,0.975)) %>%
           mutate(Size = factor(Size, levels = c("1","2","3"), labels = c("Small","Medium","Large"))))
}

# Extract psi from model (growth)
get_psi <- function(model){
  model <- as.mcmc(model)
  psi = data.frame()
  for (n in 1:2){
    for (tr in 1:4){
      for (gs in 1:2){
        for (e in 1:2){
          psi <- rbind(psi, 
                       data.frame(model[[n]][,paste("psi",as.character(gs),as.character(gs+1),"[",as.character(tr),",",as.character(e),"]", sep = "")],
                                  Size = gs,
                                  Treatment = tr,
                                  Exp_time = e))
        }
      }
    }
  }
  
  return(psi %>% 
    group_by (Size,Treatment,Exp_time) %>% 
    summarize(mean = mean(var1), min = quantile(var1,0.025), max = quantile(var1,0.975)) %>%
    mutate(Size = factor(Size, levels = c("1","2"), labels = c("Small to Medium","Medium to Large")),
           Treatment  = factor(Treatment, levels = c("1","2","3","4"), labels = Treatment_name),
           Exp_time = factor(Exp_time, levels = c("1","2"), 
                               labels = c("Before treatment","After treatment"))))
}

# Extract epsilon from model (variability of capture)
get_epsilon <- function(model){
  model <- as.mcmc(model)
  epsilon = data.frame()
  for (n in 1:2){
    for (l in 1:16){
      for (t in 1:5){
        epsilon <- rbind(epsilon, 
                         data.frame(model[[n]][,paste("epsilon","[",as.character(l),",",as.character(t),"]", sep = "")],
                                    Lake = l,
                                    Year = t))
      }
    }
  }
  
  return(epsilon %>% 
           group_by (Lake,Year) %>% 
           summarize(mean = mean(var1), min = quantile(var1,0.025), max = quantile(var1,0.975)) %>%
           mutate(Year = factor(Year, levels = c(1,2,3,4,5), labels = c(2017,2018,2019,2020,2021))))
}

# extract abundance from model
get_n <- function(model){
  model <- as.mcmc(model)
  n = data.frame()
  for (i in 1:2){
    for (tr in 1:4){
      for (gs in 1:3){
        for (t in 1:6){
          n <- rbind(n, 
                     data.frame(model[[i]][,paste("n",as.character(gs),"[",as.character(t),",",as.character(tr),"]", sep = "")],
                                Size = gs,
                                Treatment = tr,
                                Year = t))
        }
      }
    }
  }
  return (n %>% 
          group_by (Size,Treatment,Year) %>% 
          summarize(mean = mean(var1), min = quantile(var1,0.025), max = quantile(var1,0.975)) %>%
          mutate(Size = factor(Size, levels = c("1","2","3"), labels = c("Small","Medium","Large")),
                 Treatment  = factor(Treatment, levels = c("1","2","3","4"), labels = Treatment_name),
                 Year = factor(Year, levels = c(1,2,3,4,5,6), 
                                     labels = c(2016,2017,2018,2019,2020,2021))))
}

# extract deviance from model
get_deviance <- function(model){
  deviance = data.frame()
  for (n in 1:2){
    deviance <- rbind (deviance, 
                       data.frame(as.mcmc(model)[[n]][,"deviance"]))
  }
  return(deviance)
}

get_all_deviance <- function(list_model,list_name){
  deviance = data.frame()
  i = 1
  for (model in list_model){
    deviance <- rbind(deviance, data.frame(get_deviance(model),name = list_name[i]))
    i <- i+1
  }
  return(sum_deviance = deviance %>% 
           group_by(name) %>% 
           summarize(mean = mean(var1), min = quantile(var1,0.025), max = quantile(var1,0.975)))
}

# ------------------
# PLOT MODEL RESULT
# ------------------

# plot survival
plot_phi <- function(summary){
  ggplot(summary, aes(color = factor(Treatment)))+
    geom_point(aes(x = Treatment, y = mean))+
    geom_errorbar(aes(x = Treatment, ymin = min, ymax = max), width = 0.4)+
    scale_color_discrete("Treatment")+
    facet_grid(rows = vars(Size), cols = vars(Exp_time))+
    theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())+
    labs(y = "Survival probability")
}

# plot capture probability
plot_p <- function(summary){
  ggplot(summary)+
    geom_point(aes(x = Size, y = mean))+
    geom_errorbar(aes(x = Size, 
                      ymin = min, ymax = max), width = 0.4)+
    labs(y = "Capture probability", x = "Size group")+
    ylim(0,1)
}

# plot variability of capture
plot_epsilon <- function(summary){
  ggplot(summary)+
    geom_point(aes(x = factor(Lake), y = mean))+
    geom_errorbar(aes(x = factor(Lake), ymin = min, ymax = max), width = 0.4)+
    facet_wrap(~Year)+
    labs(x = "Lake", y = "Epsilon : variation around the capture probability")
}

# plot growth
plot_psi <- function(summary) {
  ggplot(summary, aes(color = factor(Treatment)))+
    geom_point(aes(x = Treatment, y = mean))+
    geom_errorbar(aes(x = Treatment, ymin = min, ymax = max), width = 5)+
    facet_grid(rows = vars(Size), cols = vars(Exp_time))+
    scale_color_discrete("Treatment")+
    theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank())+
    labs(y = "Growth probability")
}

# plot deviance
plot_deviance <- function(list_model,name_model){
  ggplot(get_all_deviance(list_model,name_model))+
    geom_point(aes(x = name, y = mean))+
    geom_errorbar(aes(x = name, ymax = max, ymin = min))+
    labs(x = "Model", y = "Deviance")+
    theme(legend.position = "none")
}

# plot abundance
plot_n <- function(summary){
  ggplot(summary)+
    geom_bar(stat = "identity", aes(x = Year, y = mean, fill = Size))+
    #geom_errorbar(aes(x = Year, ymin =min, ymax = max))+
    scale_fill_manual (values =c("#FED98E","#FE9929","#CC4C02"))+
    facet_wrap(~Treatment)
}
