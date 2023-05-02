
########## linear regression, Pearson correlation analysis, and stepwise multiple regression ##########

setwd("F:/ISAR_data_analyses-main") # Project folder location
work_dir <- getwd() # Set working directory 

library(stats) # lm function
library(tidyverse) # select function
library(Hmisc) # rcorr function
library(reshape2) # melt function
library(corrplot) # corrplot function
library(vegan) # varpart function

### Read data file 
area <- read.csv(file = "datasets/Fragment_area.csv",row.names = 1)

diversity <- read.csv(file = "diversity&soil&spatial_indices/diversity_indices.csv")
alpha <- diversity %>% select(Plant_alpha,Fungal_alpha,Bacterial_alpha)
beta <- diversity %>% select(Plant_beta,Fungal_beta,Bacterial_beta)

soil_data <- read.csv(file = "diversity&soil&spatial_indices/soil_and_spatial_indices.csv")
#soil and spatial indices were standardized to have a mean value of 0 and variance of 1#
soil_data <- soil_data[,c(3:15)]
soil_data_scale <- as.data.frame(scale(soil_data, center=T,scale=T))#数据中心化
Habitat_quality_scale <- soil_data_scale [,c(1:6)]
Heterogeneity_scale <- soil_data_scale [,c(7:13)]

rarefied_diversity<- read.csv(file = "diversity&soil&spatial_indices/rarefied_diversity.csv",row.names = 1)

all_data <- data.frame(bind_cols(area ,diversity[,-1], Habitat_quality_scale, Heterogeneity_scale, rarefied_diversity[,-1]))

########### 1.linear regression ####################################################################
##### the effects of fragment area on diversity, soil and spatial indices ##########################

diversity1 <- diversity %>% select( -Plant_Stotal_se,-Fungal_Stotal_se,-Bacterial_Stotal_se,
                                          -Plant_alpha_se,-Fungal_alpha_se,-Bacterial_alpha_se,
                                          -Plant_beta_se,-Fungal_beta_se,-Bacterial_beta_se, 
                                   -Plant_beta_within_transect_se,-Fungal_beta_within_transect_se,-Bacterial_beta_within_transect_se, )
                                          
                                           
                                 
# fragment area was log-transformed.

all_data_lm <- data.frame(bind_cols(area ,diversity1[,-1], Habitat_quality_scale, Heterogeneity_scale, rarefied_diversity[,-1]))

lm_result <- data.frame(index="", Intercept="", Slope ="",R.squared="", P.value="")
lm_result <- lm_result[-1,]

for (i in c(3:39)){
  result <- lm(all_data_lm[,i] ~ log(all_data_lm$area) )
  index <- names(all_data_lm)[i]
  Intercept<- summary(result)$coefficients[1,1]
  Slope <- summary(result)$coefficients[2,1]
  R.squared <- summary(result)$r.squared
  P.value <- summary(result)$coefficients[2,4]
  out <- data.frame(cbind(index,Intercept, Slope, R.squared, P.value))
  lm_result <- rbind(lm_result,out )
}

write.csv(lm_result, file = "results/lm_result.csv")

library(formattable)
lm_result$Intercept <- digits(lm_result$Intercept ,3)
lm_result$Slope <- digits(lm_result$Slope,3)
lm_result$R.squared <- digits(lm_result$R.squared ,2)
lm_result$P.value<- digits(lm_result$P.value ,3)
lm_result

########### 2.Pearson correlation analysis #########################################################
##### the relationships between the within-fragment (α and β) diversity of plants and microbes #####
##### with individual environmental and spatial variables ##########################################

data_within_sample <-cbind(alpha,Habitat_quality_scale)
data_among_sample <-cbind(beta,Heterogeneity_scale )

# within_sample #
within_sample <-rcorr(as.matrix(data_within_sample), type= "pearson")

W1 <- as.matrix(within_sample$r)
W1[upper.tri(W1)] <- NA
Wdf1 <- melt(as.matrix(W1), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)

W2 <- as.matrix(within_sample$P)
W2[upper.tri(W2)] <- NA
Wdf2 <- melt(as.matrix(W2), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
within_sample_pearson <-data.frame(bind_cols(Wdf1,Wdf2[,-c(1,2)]))

# among_sample #
among_sample <-rcorr(as.matrix(data_among_sample), type= "pearson")

A1 <- as.matrix(among_sample$r)
A1[upper.tri(A1)] <- NA
Adf1 <- melt(as.matrix(A1), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)

A2 <- as.matrix(among_sample$P)
A2[upper.tri(A2)] <- NA
Adf2 <- melt(as.matrix(A2), varnames = c("row", "col"), na.rm = TRUE) %>% filter(row!= col)
among_sample_pearson <-data.frame(bind_cols(Adf1,Adf2[,-c(1,2)]))

# Combine #

pearson_result <- data.frame(bind_rows(within_sample_pearson, among_sample_pearson))
colnames(pearson_result)<-c("var1","var2","r.value","P.value")

write.csv(pearson_result, "results/pearson_result.csv")

library(formattable)
pearson_result$r.value <- digits(pearson_result$r.value,2)
pearson_result$P.value <- digits(pearson_result$P.value,3)
pearson_result

########### 3.stepwise multiple regression #########################################################
###### the best model explained the within-fragment (α and β) diversity of plants and microbes #####
###### investigated whether plant diversity was an additional predictor of microbial diversity #####

### Plant_alpha ###
Plant_alpha<-lm(Plant_alpha~SM+TOC+TC+TN+pH+EC, data = all_data)
summary(step(Plant_alpha))
Plant_alpha_best_model<-lm(Plant_alpha~ TN+pH, data = all_data)
summary(Plant_alpha_best_model)

### Fungal_alpha ###
Fungal_alpha<-lm(Fungal_alpha~SM+TOC+TC+TN+pH+EC, data = all_data)
summary(step(Fungal_alpha))
Fungal_alpha_best_model<-lm(Fungal_alpha ~ TC, data = all_data)
summary(Fungal_alpha_best_model)
#adding Plant_alpha to the best model for Fungal_alpha
Fungal_alpha_additional_model<-lm(Fungal_alpha ~ Plant_alpha +TC, data = all_data)
summary(Fungal_alpha_additional_model)

### Bacterial_alpha ###
Bacterial_alpha<-lm(Bacterial_alpha~SM+TOC+TC+TN+pH+EC, data = all_data)
summary(step(Bacterial_alpha))
Bacterial_alpha_best_model<-lm(Bacterial_alpha ~ EC, data = all_data)
summary(Bacterial_alpha_best_model)
#adding Plant_alpha to the best model for Fungal_alpha
Bacterial_alpha_additional_model<-lm(Bacterial_alpha ~ Plant_alpha +EC, data = all_data)
summary(Bacterial_alpha_additional_model)

### Plant_beta ### 
Plant_beta<-lm(Plant_beta ~ Spatial_distance+SM_EU+TOC_EU+TC_EU+TN_EU+pH_EU+EC_EU, data = all_data)
summary(step(Plant_beta))
Plant_beta_best_model<-lm(Plant_beta ~ Spatial_distance, data = all_data)
summary(Plant_beta_best_model)

### Fungal_beta### 
Fungal_beta<-lm(Fungal_beta ~ Spatial_distance+SM_EU+TOC_EU+TC_EU+TN_EU+pH_EU+EC_EU, data = all_data)
summary(step(Fungal_beta))
Fungal_beta_best_model<-lm(Fungal_beta ~1, data = all_data)
summary(Fungal_beta_best_model)
#adding Plant_alpha to the best model for Fungal_alpha
Fungal_beta_additional_model<-lm(Fungal_beta ~ Plant_beta, data = all_data)
summary(Fungal_beta_additional_model)

### Bacterial_beta ###
Bacterial_beta<-lm(Bacterial_beta ~ Spatial_distance+SM_EU+TOC_EU+TC_EU+TN_EU+pH_EU+EC_EU, data = all_data)
summary(step(Bacterial_beta))
#only include pH_EU (P < 0.05)
Bacterial_beta_best_model<-lm(Bacterial_beta ~ pH_EU, data = all_data)
summary(Bacterial_beta_best_model)
#adding Plant_alpha to the best model for Fungal_alpha
Bacterial_beta_additional_model<-lm(Bacterial_beta ~ Plant_beta + pH_EU, data = all_data)
summary(Bacterial_beta_additional_model)


########### 4.variation partitioning ################################################################
###### the contribution of spatial and environmental factors to plant and microbial β diversity #####
#####################################################################################################

beta <- diversity %>% select(Plant_beta,Fungal_beta,Bacterial_beta)
Soil_Heterogeneity_scale <- soil_data_scale [,c(7:12)]
Spatial_distance_scale <- soil_data_scale [,c(13)]

corr<-cor(Soil_Heterogeneity_scale);corr

corrplot(corr,order = "original",type="upper",tl.pos = "d",tl.cex=0.6,cl.cex =0.6)
corrplot(corr,add=TRUE, type="lower", method="number",order="original",
                        number.cex = 0.6,col="black",diag=FALSE,tl.pos="n", cl.pos="n")

library(export)
graph2ppt(file="Figure/Fig_S5_Person", width=5, height=5)

# excluded heterogeneity of soil total organic C and total C, highly correlated with soil total N heterogeneity, 
# as well as conductivity heterogeneity, highly correlated with pH heterogeneity (Fig. S5). 
# Thus, heterogeneity of soil moisture, total N, and pH were used. 

Soil_Heterogeneity_scale1 <- soil_data_scale [,c(7,10,11)]

#########################
### Plant β diversity ###
#########################

Plant_vp<-varpart(beta$Plant_beta,Soil_Heterogeneity_scale1, Spatial_distance_scale)
Plant_vp
Plant_vp_plot <- plot(Plant_vp, digits = 3, Xnames = c('ENV', 'SPAT'), bg = c('#377EB8', '#E41A1C'))
Plant_vp_plot 


##########################
### Fungal β diversity ###
##########################

Fungal_vp<-varpart(beta$Fungal_beta,Soil_Heterogeneity_scale1, Spatial_distance_scale,beta$Plant_beta)
Fungal_vp
Fungal_vp_plot <-plot(Fungal_vp, digits = 3, Xnames = c('ENV', 'SPAT','PLANT'), bg = c('#377EB8', '#E41A1C','#4DAF4A'))
Fungal_vp_plot 

#############################
### Bacterial β diversity ###
#############################
Bacterial_vp<-varpart(beta$Bacterial_beta,Soil_Heterogeneity_scale1, Spatial_distance_scale,beta$Plant_beta)
Bacterial_vp
Bacterial_vp_plot <-plot(Bacterial_vp, digits = 3, Xnames = c('ENV', 'SPAT','PLANT'), bg = c('#377EB8', '#E41A1C','#4DAF4A'))
Bacterial_vp_plot 





