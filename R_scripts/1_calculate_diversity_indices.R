
########## calculated fragment-scale (Stotal), within- (α) and among-sample (β) diversity of plants, fungi and bacteria ##########

setwd("F:/ISAR_data_analyses-main") # Project folder location
work_dir <- getwd() # Set working directory 

library(iNEXT)
library(vegan)
library(tidyverse) # select function
library(plyr) # ddply function
library(reshape2) # melt function

# Read data file
Plant_table <- read.csv(file = "datasets/Plant_table.csv",row.names = 1)
Fungal_table <- read.csv(file = "datasets/Fungal_OTU_table.csv",row.names = 1)
Bacterial_table <- read.csv(file = "datasets/Bacterial_OTU_table.csv",row.names = 1)
FragmentID <- Plant_table[,1]
Fragment <- Plant_table[,2]

############ 1.calculated fragment-scale (Stotal) diversity: ##########################################
##### the estimated total plant species richness and microbial OTU richness of the whole fragment #####
#######################################################################################################

# site by species table at the gamma scale
Plant_Stotal_tab <- Plant_table[,-2]   %>% 
  group_by(FragmentID) %>%
  summarise_all(sum) %>%
  ungroup()
class(Plant_Stotal_tab) <- "data.frame"
rownames(Plant_Stotal_tab)=Plant_Stotal_tab[,1]
Plant_Stotal_tab=Plant_Stotal_tab[,-1]
Plant_Stotal_tab <- data.frame(t(Plant_Stotal_tab))

Fungal_Stotal_tab <- Fungal_table[,-2]     %>% 
  group_by(FragmentID) %>%
  summarise_all(sum) %>%
  ungroup()
class(Fungal_Stotal_tab) <- "data.frame"
rownames(Fungal_Stotal_tab)=Fungal_Stotal_tab[,1]
Fungal_Stotal_tab=Fungal_Stotal_tab[,-1]
Fungal_Stotal_tab <- data.frame(t(Fungal_Stotal_tab))

Bacterial_Stotal_tab <- Bacterial_table[,-2]    %>% 
  group_by(FragmentID) %>%
  summarise_all(sum) %>%
  ungroup()
class(Bacterial_Stotal_tab) <- "data.frame"
rownames(Bacterial_Stotal_tab)=Bacterial_Stotal_tab[,1]
Bacterial_Stotal_tab=Bacterial_Stotal_tab[,-1]
Bacterial_Stotal_tab <- data.frame(t(Bacterial_Stotal_tab))


# calculated Stotal as Chao1 estimator
set.seed(1)
Plant_S <- ChaoRichness(Plant_Stotal_tab, datatype = "abundance", conf = 0.95)
Plant_Stotal <- Plant_S %>% select(Estimator,Est_s.e.)
colnames(Plant_Stotal)<-c("Plant_Stotal","Plant_Stotal_se")

set.seed(1)
Fungal_S <- ChaoRichness(Fungal_Stotal_tab, datatype = "abundance", conf = 0.95)
Fungal_Stotal <- Fungal_S %>% select(Estimator,Est_s.e.)
colnames(Fungal_Stotal)<-c("Fungal_Stotal","Fungal_Stotal_se")

set.seed(1)
Bacterial_S <- ChaoRichness(Bacterial_Stotal_tab, datatype = "abundance", conf = 0.95)
Bacterial_Stotal <- Bacterial_S %>% select(Estimator,Est_s.e.)
colnames(Bacterial_Stotal)<-c("Bacterial_Stotal","Bacterial_Stotal_se")

Stotal <- bind_cols(Plant_Stotal,Fungal_Stotal,Bacterial_Stotal)


############## 2. Calculated within-sample (α) diversity: ###########################
##### the average Plant species richness and microbial OTU richness per quadrat #####
#####################################################################################

# Plant species richness and microbial OTU richness per quadrat
Plant_richness <- specnumber(Plant_table[,-c(1,2)])
Fungal_richness <- specnumber(Fungal_table[,-c(1,2)])
Bacterial_richness <- specnumber(Bacterial_table[,-c(1,2)])
richness_per_sample<-data.frame(cbind(Fragment,Plant_richness,Fungal_richness,Bacterial_richness))

Plant_alpha = ddply(richness_per_sample, .(Fragment), 
                    summarise,  
                    Plant_alpha =mean(Plant_richness),
                    Plant_alpha_se =sd(Plant_richness)/sqrt(length(Plant_richness)))
Fungal_alpha = ddply(richness_per_sample, .(Fragment), 
                     summarise,  
                     Fungal_alpha =mean(Fungal_richness),
                     Fungal_alpha_se =sd(Fungal_richness)/sqrt(length(Fungal_richness)))
Bacterial_alpha = ddply(richness_per_sample, .(Fragment), 
                        summarise,  
                        Bacterial_alpha =mean(Bacterial_richness),
                        Bacterial_alpha_se =sd(Bacterial_richness)/sqrt(length(Bacterial_richness)))

###### combine #####

alpha <-bind_cols(Plant_alpha[,-1],Fungal_alpha[,-1],Bacterial_alpha[,-1])



################ 3. Calculated among-sample (β) diversity: #################################################
##### the average abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats #####
############################################################################################################

############################
#### plant β diversity #####
############################

Plant_table1<- Plant_table %>% select(-FragmentID,-Fragment)
# abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats
Plant_beta1 <- data.frame(Plant_beta ="",Plant_beta_se="")
Plant_beta1 <- Plant_beta1[-1,]


for (i in 1:27) {
  if(i<=8){ tab<-Plant_table1[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Plant_table1[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Plant_table1[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Plant_table1[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Plant_table1[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  bray<- as.matrix(vegdist(tab,method = "bray"))
  bray[upper.tri(bray)] <- NA
  bray <- melt(as.matrix(bray), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- bray[,3]
  Plant_beta <- mean(a)
  Plant_beta_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Plant_beta,Plant_beta_se))
  Plant_beta1 <- rbind.data.frame(Plant_beta1,B)
  
  print(i);
  i=i+1;
}

#############################
#### Fungal β diversity #####
#############################


Fungal_table1<- Fungal_table %>% select(-FragmentID,-Fragment)
# abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats
Fungal_beta1 <- data.frame(Fungal_beta ="",Fungal_beta_se="")
Fungal_beta1 <- Fungal_beta1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-Fungal_table1[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Fungal_table1[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Fungal_table1[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Fungal_table1[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Fungal_table1[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  bray<- as.matrix(vegdist(tab,method = "bray"))
  bray[upper.tri(bray)] <- NA
  bray <- melt(as.matrix(bray), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- bray[,3]
  Fungal_beta <- mean(a)
  Fungal_beta_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Fungal_beta,Fungal_beta_se))
  Fungal_beta1 <- rbind.data.frame(Fungal_beta1,B)
  
  print(i);
  i=i+1;
}

################################
#### Bacterial β diversity #####
################################

Bacterial_table1<- Bacterial_table %>% select(-FragmentID,-Fragment)
# abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats
Bacterial_beta1 <- data.frame(Bacterial_beta ="",Bacterial_beta_se="")
Bacterial_beta1 <- Bacterial_beta1[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-Bacterial_table1[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Bacterial_table1[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Bacterial_table1[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Bacterial_table1[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Bacterial_table1[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  
  bray<- as.matrix(vegdist(tab,method = "bray"))
  bray[upper.tri(bray)] <- NA
  bray <- melt(as.matrix(bray), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  a <- bray[,3]
  Bacterial_beta <- mean(a)
  Bacterial_beta_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Bacterial_beta,Bacterial_beta_se))
  Bacterial_beta1 <- rbind.data.frame(Bacterial_beta1,B)
  
  print(i);
  i=i+1;
}


###### combine #####

beta <-bind_cols(Plant_beta1,Fungal_beta1,Bacterial_beta1)


################ 4. Calculated among-sample (β) diversity within transects: ################################
##### the average abundance-based Bray–Curtis dissimilarity among  among quadrats within transects #########
############################################################################################################

############################################
#### Plant β diversity within transect #####
############################################

Plant_table1<- Plant_table %>% select(-FragmentID,-Fragment)

# abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats
# all values of within-transect beta diversity
Plant_beta_t <- data.frame()
for (i in 1:63) {
  tab<-Plant_table1[c(i*3-2,i*3-1,i*3),]
  bray<- as.matrix(vegdist(tab,method = "bray"))
  bray[upper.tri(bray)] <- NA
  bray <- melt(as.matrix(bray), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  Plant_beta_t <- rbind.data.frame(Plant_beta_t,bray)
  print(i);
  i=i+1;
}

# the average beta within each transect within fragment (n=27)
Plant_beta_within_transect <- data.frame(Plant_beta ="",Plant_beta_se="")
Plant_beta_within_transect <- Plant_beta_within_transect[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-Plant_beta_t[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Plant_beta_t[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Plant_beta_t[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Plant_beta_t[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Plant_beta_t[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  a <- tab[,3]
  Plant_beta <- mean(a)
  Plant_beta_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Plant_beta,Plant_beta_se))
  Plant_beta_within_transect <- rbind.data.frame(Plant_beta_within_transect,B)
  
  print(i);
  i=i+1;
}

colnames(Plant_beta_within_transect) <- c('Plant_beta_within_transect','Plant_beta_within_transect_se')

#############################################
#### Fungal β diversity within transect #####
#############################################

Fungal_table1<- Fungal_table %>% select(-FragmentID,-Fragment)

# abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats
# all values of within-transect beta diversity
Fungal_beta_t <- data.frame()
for (i in 1:63) {
  tab<-Fungal_table1[c(i*3-2,i*3-1,i*3),]
  bray<- as.matrix(vegdist(tab,method = "bray"))
  bray[upper.tri(bray)] <- NA
  bray <- melt(as.matrix(bray), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  Fungal_beta_t <- rbind.data.frame(Fungal_beta_t,bray)
  print(i);
  i=i+1;
}

# the average beta within each transect within fragment (n=27)
Fungal_beta_within_transect <- data.frame(Fungal_beta ="",Fungal_beta_se="")
Fungal_beta_within_transect <- Fungal_beta_within_transect[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-Fungal_beta_t[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Fungal_beta_t[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Fungal_beta_t[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Fungal_beta_t[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Fungal_beta_t[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  a <- tab[,3]
  Fungal_beta <- mean(a)
  Fungal_beta_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Fungal_beta,Fungal_beta_se))
  Fungal_beta_within_transect <- rbind.data.frame(Fungal_beta_within_transect,B)
  
  print(i);
  i=i+1;
}

colnames(Fungal_beta_within_transect) <- c('Fungal_beta_within_transect','Fungal_beta_within_transect_se')


################################################
#### Bacterial β diversity within transect #####
################################################

Bacterial_table1<- Bacterial_table %>% select(-FragmentID,-Fragment)

# abundance-based Bray–Curtis dissimilarity among all pairwise comparisons of quadrats
# all values of within-transect beta diversity
Bacterial_beta_t <- data.frame()
for (i in 1:63) {
  tab<-Bacterial_table1[c(i*3-2,i*3-1,i*3),]
  bray<- as.matrix(vegdist(tab,method = "bray"))
  bray[upper.tri(bray)] <- NA
  bray <- melt(as.matrix(bray), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
  Bacterial_beta_t <- rbind.data.frame(Bacterial_beta_t,bray)
  print(i);
  i=i+1;
}

# the average beta within each transect within fragment (n=27)
Bacterial_beta_within_transect <- data.frame(Bacterial_beta ="",Bacterial_beta_se="")
Bacterial_beta_within_transect <- Bacterial_beta_within_transect[-1,]

for (i in 1:27) {
  if(i<=8){ tab<-Bacterial_beta_t[c(i*3-2,i*3-1,i*3),]
  
  }else if(i>8 & i<=16){
    tab<-Bacterial_beta_t[c((25+6*(i-9)):(24+6*(i-8))),]
    
  }else if(i>16 & i <= 22){
    tab<-Bacterial_beta_t[c((73+9*(i-17)):(72+9*(i-16))),]
    
  }else if(i>22 & i<=26){
    tab<-Bacterial_beta_t[c((127+12*(i-23)):(126 + 12*(i-22))),]
    
  }else if(i>26 & i<=27){
    tab<-Bacterial_beta_t[c((175+15*(i-27)):(174 + 15*(i-26))),]}
  a <- tab[,3]
  Bacterial_beta <- mean(a)
  Bacterial_beta_se <- sd(a)/sqrt(length(a))
  B <- data.frame(cbind(Bacterial_beta,Bacterial_beta_se))
  Bacterial_beta_within_transect <- rbind.data.frame(Bacterial_beta_within_transect,B)
  
  print(i);
  i=i+1;
}

colnames(Bacterial_beta_within_transect) <- c('Bacterial_beta_within_transect','Bacterial_beta_within_transect_se')

###### combine #####

beta_within_transect <-data.frame(bind_cols(Plant_beta_within_transect,
                                 Fungal_beta_within_transect,
                                 Bacterial_beta_within_transect))


#####################################################################                                                              
################### combine all diversity indices ###################
#####################################################################

diversity_indices <- data.frame(cbind (Stotal, alpha, beta, beta_within_transect))

write.csv(diversity_indices , file = "diversity&soil&spatial_indices/diversity_indices.csv")






