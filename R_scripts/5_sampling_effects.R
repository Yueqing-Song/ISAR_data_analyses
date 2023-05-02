
##################################################Sampling effects ####################################################

setwd("F:/ISAR_data_analyses-main") # Project folder location
work_dir <- getwd() # Set working directory 

library(sampling) # strata function
library(tidyverse) # select function
library(vegan)
library(reshape2) # melt function
library(plyr) # ddply function
library(mobr) # calc_biodiv function

#######################################################################################################################
###### 1. calculated diversity rarefied to an equal number of quadrats, per fragment ##################################
#######################################################################################################################

# The row and column names that will be used in the run
col_Stotal <- data.frame(matrix(1:27,nrow = 27, ncol = 1))
colnames(col_Stotal) <- c('Fragment')
col <- data.frame()
for (i in c(1:27)){
  a <- data.frame(matrix(rep(i,3),nrow = 3, ncol = 1))
  col <- rbind.data.frame(col,a)}
colnames(col) <- c('Fragment') 

##################
####   Plant   ###
##################

# The row and column names that will be used in the run
Stotal.plot <-col_Stotal
plant_Stotal_diversity <-col_Stotal

alpha.plot <-  col
plant_alpha_diversity <- col

beta.plot <-col
plant_beta_diversity <- col

# read data
ml_data <- read.csv(file = "datasets/Plant_table.csv",header = T)
L=length(unique(ml_data$FragmentID));L # number of layers
nh=rep(3,L);nh #The number of sample units at each layer

set.seed(1)
for (i in 1:1000) {  
  # srswor sampling method is used for stratified sampling of the data set, 
  # and the output result is the location of the selected sample in the population.
  sample <- strata(ml_data,stratanames = c("FragmentID"),size = nh,method="srswor");sample  
  data.strata=getdata(ml_data,sample);data.strata # Sample data is extracted
  
  ##### Stotal diversity #####
  df<- select(data.strata,-c(Sample,FragmentID,ID_unit,Prob,Stratum))
  sp.mat <- df   %>% 
    group_by(Fragment) %>%
    summarise_all(sum) %>%
    ungroup()
  class(sp.mat) <- "data.frame"
  rownames(sp.mat) = sp.mat[,1]
  sp.mat= sp.mat[,-1] 
  richness<-data.frame(specnumber(sp.mat),row.names = NULL)
  richness <-cbind(Stotal.plot,richness)
  names(richness)<-c("Fragment", "S") 
  plant_Stotal_diversity  <-left_join(plant_Stotal_diversity  ,richness,by="Fragment")
  
  ##### alpha diversity #####
  df1<- select(data.strata,-c(Fragment,FragmentID,ID_unit,Prob,Stratum))
  row.names(df1) = NULL
  rownames(df1)=df1[,1]
  df1=df1[,-1]
  head(df1)
  A<-data.frame(specnumber(df1),row.names = NULL)
  names(A)<-c("") 
  plant_alpha_diversity<-cbind(plant_alpha_diversity,A)
  
  ##### diversity ####
  beta <- data.frame()
  for (j in 1:27) {
    B<-df1[c(j*3-2,j*3-1,j*3),]
    b1 <- vegdist(B,method = "bray")
    b1 <- as.matrix(b1)
    b1[upper.tri(b1)] <- NA
    BB <- melt(as.matrix(b1), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
    beta <-rbind(beta,BB);
    j=j+1
  }
  BBB<- select(beta,-c(row,col))
  names(BBB)<-c( "") 
  plant_beta_diversity<-cbind(plant_beta_diversity,BBB)
  
  print(i);
  i=i+1;
}

# Calculate rarefied_diversity
# the average of the diversity results across a thousand samples

rarefied_plant_Stotal<- data.frame(rowMeans(plant_Stotal_diversity[,-1]))
colnames(rarefied_plant_Stotal)<- 'rarefied_plant_Stotal'

d=ddply(plant_alpha_diversity, .(Fragment), function(x) colMeans(x[-1]))
rarefied_plant_alpha<- data.frame(rowMeans(d[,-1]))
colnames(rarefied_plant_alpha)<- 'rarefied_plant_alpha'

f=ddply( plant_beta_diversity, .(Fragment), function(x) colMeans(x[-1]))
rarefied_plant_beta<- data.frame(rowMeans(f[,-1]))
colnames(rarefied_plant_beta) <- 'rarefied_plant_beta'

rarefied_plant_diversity_plot <- cbind(rarefied_plant_Stotal,rarefied_plant_alpha,rarefied_plant_beta)


##################
####   Fungi   ###
##################

# The row and column names that will be used in the run
Stotal.plot <-col_Stotal
fungal_Stotal_diversity <-col_Stotal

alpha.plot <-  col
fungal_alpha_diversity <- col

beta.plot <-col
fungal_beta_diversity <- col

# read data
ml_data <- read.csv(file = "datasets/Fungal_OTU_table.csv",header = T)
L=length(unique(ml_data$FragmentID));L # number of layers
nh=rep(3,L);nh #The number of sample units at each layer

set.seed(1)
for (i in 1:1000) {  
  # srswor sampling method is used for stratified sampling of the data set, 
  # and the output result is the location of the selected sample in the population.
  sample <- strata(ml_data,stratanames = c("FragmentID"),size = nh,method="srswor");sample  
  data.strata=getdata(ml_data,sample);data.strata # Sample data is extracted
  
  ##### Stotal diversity #####
  df<- select(data.strata,-c(Sample,FragmentID,ID_unit,Prob,Stratum))
  sp.mat <- df   %>% 
    group_by(Fragment) %>%
    summarise_all(sum) %>%
    ungroup()
  class(sp.mat) <- "data.frame"
  rownames(sp.mat) = sp.mat[,1]
  sp.mat= sp.mat[,-1] 
  richness<-data.frame(specnumber(sp.mat),row.names = NULL)
  richness <-cbind(Stotal.plot,richness)
  names(richness)<-c("Fragment", "S") 
  fungal_Stotal_diversity  <-left_join(fungal_Stotal_diversity  ,richness,by="Fragment")
  
  ##### alpha diversity #####
  df1<- select(data.strata,-c(Fragment,FragmentID,ID_unit,Prob,Stratum))
  row.names(df1) = NULL
  rownames(df1)=df1[,1]
  df1=df1[,-1]
  head(df1)
  A<-data.frame(specnumber(df1),row.names = NULL)
  names(A)<-c("") 
  fungal_alpha_diversity<-cbind(fungal_alpha_diversity,A)
  
  ##### diversity ####
  beta <- data.frame()
  for (j in 1:27) {
    B<-df1[c(j*3-2,j*3-1,j*3),]
    b1 <- vegdist(B,method = "bray")
    b1 <- as.matrix(b1)
    b1[upper.tri(b1)] <- NA
    BB <- melt(as.matrix(b1), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
    beta <-rbind(beta,BB);
    j=j+1
  }
  BBB<- select(beta,-c(row,col))
  names(BBB)<-c( "") 
  fungal_beta_diversity<-cbind(fungal_beta_diversity,BBB)
  
  print(i);
  i=i+1;
}

# Calculate rarefied_diversity
# the average of the diversity results across a thousand samples

rarefied_fungal_Stotal<- data.frame(rowMeans(fungal_Stotal_diversity[,-1]))
colnames(rarefied_fungal_Stotal)<- 'rarefied_fungal_Stotal'

d=ddply(fungal_alpha_diversity, .(Fragment), function(x) colMeans(x[-1]))
rarefied_fungal_alpha<- data.frame(rowMeans(d[,-1]))
colnames(rarefied_fungal_alpha)<- 'rarefied_fungal_alpha'

f=ddply( fungal_beta_diversity, .(Fragment), function(x) colMeans(x[-1]))
rarefied_fungal_beta<- data.frame(rowMeans(f[,-1]))
colnames(rarefied_fungal_beta) <- 'rarefied_fungal_beta'

rarefied_fungal_diversity_plot <- cbind(rarefied_fungal_Stotal,rarefied_fungal_alpha,rarefied_fungal_beta)

###################
####  bacteria  ###
###################

# The row and column names that will be used in the run
Stotal.plot <-col_Stotal
bacterial_Stotal_diversity <-col_Stotal

alpha.plot <-  col
bacterial_alpha_diversity <- col

beta.plot <-col
bacterial_beta_diversity <- col

# read data
ml_data <- read.csv(file = "datasets/Bacterial_OTU_table.csv",header = T)
L=length(unique(ml_data$FragmentID));L # number of layers
nh=rep(3,L);nh #The number of sample units at each layer

set.seed(1)
for (i in 1:1000) {  
  # srswor sampling method is used for stratified sampling of the data set, 
  # and the output result is the location of the selected sample in the population.
  sample <- strata(ml_data,stratanames = c("FragmentID"),size = nh,method="srswor");sample  
  data.strata=getdata(ml_data,sample);data.strata # Sample data is extracted
  
  ##### Stotal diversity #####
  df<- select(data.strata,-c(Sample,FragmentID,ID_unit,Prob,Stratum))
  sp.mat <- df   %>% 
    group_by(Fragment) %>%
    summarise_all(sum) %>%
    ungroup()
  class(sp.mat) <- "data.frame"
  rownames(sp.mat) = sp.mat[,1]
  sp.mat= sp.mat[,-1] 
  richness<-data.frame(specnumber(sp.mat),row.names = NULL)
  richness <-cbind(Stotal.plot,richness)
  names(richness)<-c("Fragment", "S") 
  bacterial_Stotal_diversity  <-left_join(bacterial_Stotal_diversity  ,richness,by="Fragment")
  
  ##### alpha diversity #####
  df1<- select(data.strata,-c(Fragment,FragmentID,ID_unit,Prob,Stratum))
  row.names(df1) = NULL
  rownames(df1)=df1[,1]
  df1=df1[,-1]
  head(df1)
  A<-data.frame(specnumber(df1),row.names = NULL)
  names(A)<-c("") 
  bacterial_alpha_diversity<-cbind(bacterial_alpha_diversity,A)
  
  ##### diversity ####
  beta <- data.frame()
  for (j in 1:27) {
    B<-df1[c(j*3-2,j*3-1,j*3),]
    b1 <- vegdist(B,method = "bray")
    b1 <- as.matrix(b1)
    b1[upper.tri(b1)] <- NA
    BB <- melt(as.matrix(b1), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
    beta <-rbind(beta,BB);
    j=j+1
  }
  BBB<- select(beta,-c(row,col))
  names(BBB)<-c( "") 
  bacterial_beta_diversity<-cbind(bacterial_beta_diversity,BBB)
  
  print(i);
  i=i+1;
}

# Calculate rarefied_diversity
# the average of the diversity results across a thousand samples

rarefied_bacterial_Stotal<- data.frame(rowMeans(bacterial_Stotal_diversity[,-1]))
colnames(rarefied_bacterial_Stotal)<- 'rarefied_bacterial_Stotal'

d=ddply(bacterial_alpha_diversity, .(Fragment), function(x) colMeans(x[-1]))
rarefied_bacterial_alpha<- data.frame(rowMeans(d[,-1]))
colnames(rarefied_bacterial_alpha)<- 'rarefied_bacterial_alpha'

f=ddply( bacterial_beta_diversity, .(Fragment), function(x) colMeans(x[-1]))
rarefied_bacterial_beta<- data.frame(rowMeans(f[,-1]))
colnames(rarefied_bacterial_beta) <- 'rarefied_bacterial_beta'

rarefied_bacterial_diversity_plot <- cbind(rarefied_bacterial_Stotal,rarefied_bacterial_alpha,rarefied_bacterial_beta)


###################
####  combine   ###
###################

rarefied_diversity_plot <- cbind(rarefied_plant_diversity_plot,
                                 rarefied_fungal_diversity_plot,
                                 rarefied_bacterial_diversity_plot)

##########################################################################################################################
######## 2. calculated diversity rarefied to an equal number of individuals/sequences, per fragment ######################
##########################################################################################################################

##################
####   Plant   ###
##################

dat<-read.csv(file = "datasets/Plant_table.csv",header = T)
dat <- dat %>%  select(-FragmentID, -Sample)

# site by species table at the gamma scale
Stotal_tab <- dat  %>% 
  group_by(Fragment) %>%
  summarise_all(sum) %>%
  ungroup()
class(Stotal_tab) <- "data.frame"

# rarefied to minimum number of individuals/sequences 
n_sites <- rowSums(Stotal_tab[,-1])
n_ref <- min(n_sites)
set.seed(1)
plant_Srare <- calc_biodiv(Stotal_tab[,-1], 
                           groups = Stotal_tab$Fragment,
                           index = c("S_n"),
                           effort = n_ref,
                           extrapolate = T,
                           return_NA = F)
colnames(plant_Srare)[4] <- 'plant_Srare'
plant_Srare <- plant_Srare %>%  select( -index, -effort)

##################
####   Fungi   ###
##################

dat<-read.csv("datasets/Fungal_OTU_table.csv",header = T)
dat <- dat %>%  select(-FragmentID, -Sample)

# site by species table at the gamma scale
Stotal_tab <- dat  %>% 
  group_by(Fragment) %>%
  summarise_all(sum) %>%
  ungroup()
class(Stotal_tab) <- "data.frame"

# rarefied to minimum number of individuals/sequences 
n_sites <- rowSums(Stotal_tab[,-1])
n_ref <- min(n_sites)
set.seed(1)
fungal_Srare <- calc_biodiv(Stotal_tab[,-1], 
                            groups = Stotal_tab$Fragment,
                            index = c("S_n"),
                            effort = n_ref,
                            extrapolate = T,
                            return_NA = F)
colnames(fungal_Srare)[4] <- 'fungal_Srare'
fungal_Srare <- fungal_Srare %>%  select(-group, -index, -effort)

###################
####  bacteria  ###
###################

dat<-read.csv("datasets/Bacterial_OTU_table.csv",header = T)
dat <- dat %>% select(-FragmentID, -Sample)

# site by species table at the gamma scale
Stotal_tab <- dat  %>% 
  group_by(Fragment) %>%
  summarise_all(sum) %>%
  ungroup()
class(Stotal_tab) <- "data.frame"

# rarefied to minimum number of individuals/sequences 
n_sites <- rowSums(Stotal_tab[,-1])
n_ref <- min(n_sites)
set.seed(1)
bacterial_Srare <- calc_biodiv(Stotal_tab[,-1], 
                               groups = Stotal_tab$Fragment,
                               index = c("S_n"),
                               effort = n_ref,
                               extrapolate = T,
                               return_NA = F)
colnames(bacterial_Srare)[4] <- 'bacterial_Srare'
bacterial_Srare <- bacterial_Srare %>%  select(-group, -index, -effort)

Srare <- cbind(plant_Srare,fungal_Srare,bacterial_Srare)


rarefied_diversity <- cbind(Srare,rarefied_diversity_plot) 
colnames(rarefied_diversity)[1] <- 'Fragment'

write.csv(rarefied_diversity, file = "diversity&soil&spatial_indices/rarefied_diversity.csv")

