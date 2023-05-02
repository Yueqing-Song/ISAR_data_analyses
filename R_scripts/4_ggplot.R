
################################################ make figures #######################################################

setwd("F:/ISAR_data_analyses") # Project folder location
work_dir <- getwd() # Set working directory 

library(tidyverse) 
library(ggplot2)
library(cowplot)

#install.packages("Cairo")
area <- read.csv(file = "datasets/Fragment_area.csv",row.names = 1)
diversity <- read.csv(file = "diversity&soil&spatial_indices/diversity_indices.csv")
rarefied_diversity<- read.csv(file = "diversity&soil&spatial_indices/rarefied_diversity.csv",row.names = 1)

data <- data.frame(bind_cols(area ,diversity[,-1], rarefied_diversity[,-1]))

############# fragment-scale (Stotal), within- (α) and among-sample (β) diversity #########################

### 1.Plant_Stotal ~log(area) ###
theme_set(theme_bw())
Plant_Stotal <- ggplot(data, mapping = aes(x =log(area), y = Plant_Stotal))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_errorbar(aes(ymin=Plant_Stotal-Plant_Stotal_se, ymax=Plant_Stotal+Plant_Stotal_se), 
                width=0.15, linewidth = 0.2,colour = "#4DAF4A")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(bold(paste("log Area (m" ^2,")"))), y="Stotal\n(Estimated richness per fragment)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Plant_Stotal

### 2.Fungal_Stotal ~log(area) ###
Fungal_Stotal <- ggplot(data, mapping = aes(x =log(area), y = Fungal_Stotal))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_errorbar(aes(ymin=Fungal_Stotal-Fungal_Stotal_se, ymax=Fungal_Stotal+Fungal_Stotal_se), 
                width=0.15, linewidth  = 0.2,colour = "#377EB8")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth  = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(bold(paste("log Area (m" ^2,")"))), y="Stotal\n(Estimated richness per fragment)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =0.4, linetype="solid"))
Fungal_Stotal

### 3.Bacterial_Stotal ~log(area) ###
Bacterial_Stotal <- ggplot(data, mapping = aes(x =log(area), y = Bacterial_Stotal))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  geom_errorbar(aes(ymin=Bacterial_Stotal-Bacterial_Stotal_se, ymax=Bacterial_Stotal+Bacterial_Stotal_se), 
                width=0.15, linewidth = 0.2,colour = "#FF7F00")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(bold(paste("log Area (m" ^2,")"))), y="Stotal\n(Estimated richness per fragment)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Bacterial_Stotal

### 4.Plant_Stotal ~log(area) ###

Plant_alpha <- ggplot(data, mapping = aes(x =log(area), y = Plant_alpha))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_errorbar(aes(ymin=Plant_alpha-Plant_alpha_se, ymax=Plant_alpha+Plant_alpha_se), 
                width=0.15, linewidth = 0.2,colour = "#4DAF4A")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth= 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Alpha diversity\n(Richness per sample)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Plant_alpha

### 5.Fungal_alpha ~log(area) ###
Fungal_alpha <- ggplot(data, mapping = aes(x =log(area), y = Fungal_alpha))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_errorbar(aes(ymin=Fungal_alpha-Fungal_alpha_se, ymax=Fungal_alpha+Fungal_alpha_se), 
                width=0.15, linewidth = 0.2,colour = "#377EB8")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Alpha diversity\n(Richness per sample)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Fungal_alpha

  
### 6.Bacterial_alpha ~log(area) ###
Bacterial_alpha <- ggplot(data, mapping = aes(x =log(area), y = Bacterial_alpha))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  geom_errorbar(aes(ymin=Bacterial_alpha-Bacterial_alpha_se, ymax=Bacterial_alpha+Bacterial_alpha_se), 
                width=0.15, linewidth = 0.2,colour = "#FF7F00")+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Alpha diversity\n(Richness per sample)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Bacterial_alpha

### 7.Plant_beta ~log(area) ###
Plant_beta <- ggplot(data, mapping = aes(x =log(area), y = Plant_beta))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_errorbar(aes(ymin=Plant_beta-Plant_beta_se, ymax=Plant_beta+Plant_beta_se), 
                width=0.15,linewidth = 0.2,colour = "#4DAF4A")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Beta diversity\n(Bray-Curtis dissimilarity)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Plant_beta

### 8.Fungal_beta ~log(area) ###
Fungal_beta <- ggplot(data, mapping = aes(x =log(area), y = Fungal_beta))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_errorbar(aes(ymin=Fungal_beta-Fungal_beta_se, ymax=Fungal_beta+Fungal_beta_se), 
                width=0.15, linewidth = 0.2,colour = "#377EB8")+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Beta diversity\n(Bray-Curtis dissimilarity)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Fungal_beta

### 9.Bacterial_beta ~log(area) ###
Bacterial_beta <- ggplot(data, mapping = aes(x =log(area), y = Bacterial_beta))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  geom_errorbar(aes(ymin=Bacterial_beta-Bacterial_beta_se, ymax=Bacterial_beta+Bacterial_beta_se), 
                width=0.15, linewidth = 0.2,colour = "#FF7F00")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Beta diversity\n(Bray-Curtis dissimilarity)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Bacterial_beta


### combine ###

Fig_2_diversity <- plot_grid(Plant_Stotal,Fungal_Stotal,Bacterial_Stotal, 
                           Plant_alpha,Fungal_alpha,Bacterial_alpha,
                           Plant_beta,Fungal_beta,Bacterial_beta,
                           labels=c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                           label_size=8, label_fontface = "bold", 
                           hjust=-6.5, vjust=1, ncol=3, align="hv")
Fig_2_diversity

ggsave(filename = "Figure/Fig_2_diversity.png", plot = Fig_2_diversity, 
       device = "png",dpi=500, width = 20, height = 20, units = "cm")


library(export)
graph2ppt(file="Figure/Fig_2_diversity", width=7, height=6)


######################## diversity rarefied to an equal number of individuals/sequences ####################

### 1. plant_Srare  ~log(area) ###
plant_Srare <- ggplot(data, mapping = aes(x =log(area), y = plant_Srare ))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth= 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied Stotal") +
  theme(aspect.ratio = 0.8,
        axis.title.x =  element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
plant_Srare

### 2.fungal_Srare ~log(area) ###
fungal_Srare <- ggplot(data, mapping = aes(x =log(area), y = fungal_Srare))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied Stotal") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
fungal_Srare

### 3.bacterial_Srare ~log(area) ###
bacterial_Srare <- ggplot(data, mapping = aes(x =log(area), y = bacterial_Srare))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied Stotal") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
bacterial_Srare


### combine ###

Fig_S2_diversity <- plot_grid(plant_Srare,fungal_Srare,bacterial_Srare, 
                              labels=c("a", "b", "c"),
                              label_size=8, label_fontface = "bold", 
                              hjust=-6.5, vjust=1, ncol=3, align="hv")
Fig_S2_diversity

ggsave(filename = "Figure/Fig_S2_diversity.png", plot = Fig_S2_diversity, 
       device = "png",dpi=500, width = 20, height = 7, units = "cm")


library(export)
graph2ppt(file="Figure/Fig_S2_diversity", width=7, height=3)


#################### diversity rarefied to an equal number of individuals/sequences, per fragment  #################### 

### 1.rarefied_plant_Stotal ~log(area) ###
rarefied_plant_Stotal <- ggplot(data, mapping = aes(x =log(area), y = rarefied_plant_Stotal))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(bold(paste("log Area (m" ^2,")"))), y="Rarefied Stotal") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_plant_Stotal

### 2.rarefied_fungal_Stotal ~log(area) ###
rarefied_fungal_Stotal <- ggplot(data, mapping = aes(x =log(area), y = rarefied_fungal_Stotal))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth  = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(bold(paste("log Area (m" ^2,")"))), y="Rarefied Stotal)") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =0.4, linetype="solid"))
rarefied_fungal_Stotal

### 3.rarefied_bacterial_Stotal ~log(area) ###
rarefied_bacterial_Stotal <- ggplot(data, mapping = aes(x =log(area), y = rarefied_bacterial_Stotal))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(bold(paste("log Area (m" ^2,")"))), y="Rarefied Stotal") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_bacterial_Stotal

### 4.rarefied_plant_alpha ~log(area) ###

rarefied_plant_alpha <- ggplot(data, mapping = aes(x =log(area), y = rarefied_plant_alpha))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth= 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied alpha diversity") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_plant_alpha

### 5.rarefied_fungal_alpha ~log(area) ###
rarefied_fungal_alpha <- ggplot(data, mapping = aes(x =log(area), y = rarefied_fungal_alpha))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied alpha diversity") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_fungal_alpha


### 6.rarefied_bacterial_alpha ~log(area) ###
rarefied_bacterial_alpha <- ggplot(data, mapping = aes(x =log(area), y = rarefied_bacterial_alpha))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied alpha diversity") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_bacterial_alpha

### 7.rarefied_plant_beta ~log(area) ###
rarefied_plant_beta <- ggplot(data, mapping = aes(x =log(area), y = rarefied_plant_beta))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied beta diversity") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_plant_beta

### 8.rarefied_fungal_beta ~log(area) ###
rarefied_fungal_beta <- ggplot(data, mapping = aes(x =log(area), y = rarefied_fungal_beta))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied beta diversity") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_fungal_beta

### 9.rarefied_bacterial_beta ~log(area) ###
rarefied_bacterial_beta <- ggplot(data, mapping = aes(x =log(area), y = rarefied_bacterial_beta))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Rarefied beta diversity") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
rarefied_bacterial_beta


### combine ###

Fig_S3_diversity <- plot_grid(rarefied_plant_Stotal,rarefied_fungal_Stotal,rarefied_bacterial_Stotal, 
                             rarefied_plant_alpha,rarefied_fungal_alpha,rarefied_bacterial_alpha,
                             rarefied_plant_beta,rarefied_fungal_beta,rarefied_bacterial_beta,
                             labels=c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                             label_size=8, label_fontface = "bold", 
                             hjust=-6.5, vjust=1, ncol=3, align="hv")
Fig_S3_diversity

ggsave(filename = "Figure/Fig_S3_diversity.png", plot = Fig_S3_diversity, 
       device = "png",dpi=500, width = 20, height = 20, units = "cm")


library(export)
graph2ppt(file="Figure/Fig_S3_diversity", width=7, height=6)



############################# among-sample (β) diversity within transects ########################################

theme_set(theme_bw())

### 1.Plant_beta_within_transect ~log(area) ###
Plant_beta_within_transect <- ggplot(data, mapping = aes(x =log(area), y = Plant_beta_within_transect))+
  geom_point(color = "#4DAF4A",size=1,shape=19)+
  geom_errorbar(aes(ymin=Plant_beta_within_transect-Plant_beta_within_transect_se, 
                    ymax=Plant_beta_within_transect+Plant_beta_within_transect_se), 
                width=0.15, linewidth = 0.2,colour = "#4DAF4A")+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Beta diversity within transect") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Plant_beta_within_transect

### 2.Fungal_beta_within_transect ~log(area) ###
Fungal_beta_within_transect <- ggplot(data, mapping = aes(x =log(area), y = Fungal_beta_within_transect))+
  geom_point(color = "#377EB8",size=1,shape=19)+
  geom_errorbar(aes(ymin=Fungal_beta_within_transect-Fungal_beta_within_transect_se, 
                    ymax=Fungal_beta_within_transect+Fungal_beta_within_transect_se), 
                width=0.15, linewidth = 0.2,colour = "#377EB8")+
  #geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Beta diversity within transect") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Fungal_beta_within_transect

### 3.Bacterial_beta_within_transect ~log(area) ###
Bacterial_beta_within_transect <- ggplot(data, mapping = aes(x =log(area), y = Bacterial_beta_within_transect))+
  geom_point(color = "#FF7F00",size=1,shape=19)+
  geom_errorbar(aes(ymin=Bacterial_beta_within_transect-Bacterial_beta_within_transect_se, 
                    ymax=Bacterial_beta_within_transect+Bacterial_beta_within_transect_se), 
                width=0.15, linewidth = 0.2,colour = "#FF7F00")+
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE,linewidth = 0.5,colour = "black")+
  coord_cartesian(xlim=c(1.9,8.2)) +  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x=expression(paste("log Area (m" ^2,")")), y="Beta diversity within transect") +
  theme(aspect.ratio = 0.8,
        axis.title.x = element_text(size = 8,hjust=0.5,vjust=0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 6),
        axis.text.y = element_text(colour = "black",size = 6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.4, linetype="solid"))
Bacterial_beta_within_transect


### combine ###

Fig_S4_beta_within_transect <- plot_grid(Plant_beta_within_transect ,
                                         Fungal_beta_within_transect ,
                                         Bacterial_beta_within_transect ,
                             labels=c("a", "b", "c"),
                             label_size=8, label_fontface = "bold", 
                             hjust=-6.5, vjust=1, ncol=3, align="hv")
Fig_S4_beta_within_transect

ggsave(filename = "Figure/Fig_S4_beta_within_transect.png", plot = Fig_S4_beta_within_transect, 
       device = "png",dpi=500, width = 20, height = 7, units = "cm")


library(export)
graph2ppt(file="Figure/Fig_S4_beta_within_transect", width=7, height=3)

