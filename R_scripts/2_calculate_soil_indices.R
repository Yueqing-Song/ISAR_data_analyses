
########## calculated habitat quality, habitat heterogeneity, and spatial distance on each fragment ##########

setwd("F:/ISAR_data_analyses") # Project folder location
work_dir <- getwd() # Set working directory 

library(vegan)
library(tidyverse) # select function
library(plyr) # ddply function
library(reshape2) # melt function
#install.packages("devtools")
#devtools::install_github("renkun-ken/formattable")
library(formattable) # digits function , Keep the last four digits of the decimal

data <-  read.csv(file = "datasets/Soil_properties&spatial.csv")

############### 1. calculate habitat quality on each fragment: ############################
##### the average individual soil properties per quadrat within a fragment ################
###########################################################################################

soil_data <- data %>% select(-Sample, - Transect, -Quadrat,-X,-Y)
Habitat_quality <- ddply(soil_data, .(Fragment), function(x) colMeans(x[-1]))


############## 2. calculate habitat heterogeneity and spatial distance on each fragment: ############################################
##### the average dissimilarity in individual soil properties (based on Euclidean distances) between quadrats within a fragment #####
##### the average spatial distance among quadrats within each fragment ##############################################################


########################################
### soil moisture (SM) heterogeneity ###
########################################

SM_data <- data %>%  select(Sample,Fragment,SM)
rownames(SM_data) = SM_data[,1]
SM_data = SM_data[,-1] 
SM_data <- SM_data %>%  select(-Fragment)

SM_EU1 <- data.frame(SM_EU ="",SM_EU_se="")
SM_EU1 <- SM_EU1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-SM_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-SM_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-SM_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-SM_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-SM_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  SM_EU <- mean(a)
  SM_EU_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(SM_EU,SM_EU_se))
  SM_EU1 <- rbind.data.frame(SM_EU1,B)
  
  print(i);
  i=i+1;
}

################################################
### total organic carbon (TOC) heterogeneity ###
################################################

TOC_data <- data %>% select(Sample,Fragment,TOC)
rownames(TOC_data) = TOC_data[,1]
TOC_data = TOC_data[,-1] 
TOC_data <- TOC_data %>% select(-Fragment)

TOC_EU1 <- data.frame(TOC_EU ="",TOC_EU_se="")
TOC_EU1 <- TOC_EU1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-TOC_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-TOC_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-TOC_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-TOC_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-TOC_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  TOC_EU <- mean(a)
  TOC_EU_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(TOC_EU,TOC_EU_se))
  TOC_EU1 <- rbind.data.frame(TOC_EU1,B)
  
  print(i);
  i=i+1;
}

#######################################
### total carbon (TC) heterogeneity ###
#######################################

TC_data <- data %>%  select(Sample,Fragment,TC)
rownames(TC_data) = TC_data[,1]
TC_data = TC_data[,-1]
TC_data <- TC_data %>% select(-Fragment)

TC_EU1 <- data.frame(TC_EU ="",TC_EU_se="")
TC_EU1 <- TC_EU1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-TC_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-TC_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-TC_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-TC_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-TC_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  TC_EU <- mean(a)
  TC_EU_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(TC_EU,TC_EU_se))
  TC_EU1 <- rbind.data.frame(TC_EU1,B)
  
  print(i);
  i=i+1;
}

#########################################
### total nitrogen (TN) heterogeneity ###
#########################################

# total organic carbon (TOC) euclidean distance among all pairwise comparisons of quadrats
TN_data <- data %>%  select(Sample,Fragment,TN)
rownames(TN_data) = TN_data[,1]
TN_data = TN_data[,-1] 
TN_data <- TN_data %>%  select(-Fragment)

TN_EU1 <- data.frame(TN_EU ="",TN_EU_se="")
TN_EU1 <- TN_EU1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-TN_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-TN_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-TN_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-TN_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-TN_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  TN_EU <- mean(a)
  TN_EU_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(TN_EU,TN_EU_se))
  TN_EU1 <- rbind.data.frame(TN_EU1,B)
  
  print(i);
  i=i+1;
}

########################
### pH heterogeneity ###
########################

# pH euclidean distance among all pairwise comparisons of quadrats
pH_data <- data %>% select(Sample,Fragment,pH)
rownames(pH_data) = pH_data[,1]
pH_data = pH_data[,-1] 
pH_data <- pH_data %>% select(-Fragment)

pH_EU1 <- data.frame(pH_EU ="",pH_EU_se="")
pH_EU1 <- pH_EU1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-pH_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-pH_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-pH_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-pH_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-pH_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  pH_EU <- mean(a)
  pH_EU_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(pH_EU,pH_EU_se))
  pH_EU1 <- rbind.data.frame(pH_EU1,B)
  
  print(i);
  i=i+1;
}

##################################################
### electrical conductivity (EC) heterogeneity ###
##################################################

EC_data <- data %>% select(Sample,Fragment,EC)
rownames(EC_data) = EC_data[,1]
EC_data = EC_data[,-1] 
EC_data <- EC_data %>% select(-Fragment)

EC_EU1 <- data.frame(EC_EU ="",EC_EU_se="")
EC_EU1 <- EC_EU1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-EC_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-EC_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-EC_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-EC_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-EC_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  EC_EU <- mean(a)
  EC_EU_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(EC_EU,EC_EU_se))
  EC_EU1 <- rbind.data.frame(EC_EU1,B)
  
  print(i);
  i=i+1;
}

########################
### Spatial distance ###
########################
# Spatial distance-euclidean distance among all pairwise comparisons of quadrats
# x and Y were coordinates of samples in the Projection Coordinate System, Unit: meter.

Spatial_data <- data %>% select(Sample,Fragment,X,Y)
rownames(Spatial_data) = Spatial_data[,1]
Spatial_data = Spatial_data[,-1] 
Spatial_data <- Spatial_data %>% select(-Fragment)

#Keep the last four digits of the decimal
Spatial_data$X <- digits(Spatial_data$X ,4)
Spatial_data$Y <- digits(Spatial_data$Y ,4)
Spatial_data <- data.frame(Spatial_data)

Spatial_distance1 <- data.frame(Spatial_distance ="",Spatial_distance_se="")
Spatial_distance1 <- Spatial_distance1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-Spatial_data[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Spatial_data[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Spatial_data[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Spatial_data[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Spatial_data[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  EU<- as.matrix(vegdist(tab,method = "euclidean"))
  EU[upper.tri(EU)] <- NA
  EU <- melt(as.matrix(EU), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- EU[,3]
  Spatial_distance <- mean(a)
  Spatial_distance_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Spatial_distance,Spatial_distance_se))
  Spatial_distance1 <- rbind.data.frame(Spatial_distance1,B)
  
  print(i);
  i=i+1;
}

###### combine #####

soil_EU <-data.frame(bind_cols(SM_EU1, TOC_EU1,TC_EU1,TN_EU1,pH_EU1,EC_EU1,Spatial_distance1))

soil_EU <- soil_EU %>% select(SM_EU, TOC_EU,TC_EU,TN_EU,pH_EU,EC_EU,Spatial_distance)


#############################################                                                            
### combine all soil and spatial indices ####
#############################################

soil_and_spatial_indices <- data.frame(cbind (Habitat_quality, soil_EU))

write.csv(soil_and_spatial_indices, file = "diversity&soil&spatial_indices/soil_and_spatial_indices.csv")


