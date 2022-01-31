### Nutrient Regime Clustering Analysis
# Created by Jamie Kerlin 
# Created on 2022_01_30
##########################################################################

### Load libraries #######################################################
library(tidyverse)
library(here)
library(VarSelLCM)
library(ggmap)
library(PNWColors)

register_google(key = "")

### Load data ############################################################
data <- read_csv("https://raw.githubusercontent.com/njsilbiger/NutrientRegimes/main/Data/NutrientAll.csv")
data_na <- read_csv("https://raw.githubusercontent.com/njsilbiger/NutrientRegimes/main/Data/NutrientAll_na.csv")
metadata <- read_csv(here("Data/islandwide_sites_allMetadata.csv"))  


### Remove large silicate points #########################################
data_na1 <- data_na %>% 
  filter(Silicate_May < 1.5,
         Silicate_August < 1.5)
             
### Join data and metadata ###############################################
metadata <- metadata %>% #rename site column for matching columns
  rename("Site_Number" = Site) %>%
  select(c(Site_Number, Island_shore, Habitat)) #only select needed columns

full_data <- left_join(data_na1, metadata) %>% #join data and metadata
  unite("Shore_Habitat", Island_shore:Habitat, sep = "_", na.rm = TRUE) #unite shore and habitat

#new object of class of known characteristic
ztrue <- full_data$Shore_Habitat #here, using shore and habitat as known characteristic
view(ztrue)

#remove columns not needed
x <- full_data %>% 
  select(!c(Site_Number, Lat, Lon, Date_May, Time_May, Shore_Habitat))

xlog1 <- x %>% #log transform 
  mutate_all(.funs = ~log(.x+0.1)) #log transform

xlog_scale <- scale(xlog1, center = TRUE, scale = TRUE) #scale data


#Cluster without variable selection
res_without <- VarSelCluster(xlog_scale, gvals = 1:8, nbcores = 4, vbleSelec = FALSE, crit.varsel = "BIC")
BIC(res_without) #Bayesian Information Criteria
print(res_without) #more detailed output

#Cluster with variable selection
res_with <- VarSelCluster(xlog_scale, gvals = 1:8, nbcores = 4, vbleSelec = TRUE, crit.varsel = "BIC")
BIC(res_with)
print(res_with) #more detailed output

#Compare partition accuracy using ARI- computed between true partition (ztrue) and its estimators
# 0 = partitions are independent, 1 = partitions are equals
# variable selection improves ARI- ARI not used for model selection
ARI(ztrue, fitted(res_without))
ARI(ztrue, fitted(res_with))

#Estimated partition
fitted(res_with)

#Estimated probabilities of classification
head(fitted(res_with, type = "probability"))

#Summary
summary(res_with)

#Discriminative power of the variables
plot(res_with)


### Save the clusters as a column in new object

full_data1 <- full_data
full_data1$cluster <- fitted(res_with)

full_data1 <- full_data1 %>%
  mutate("cluster" = factor(cluster))

### Map clusters 
moorea <- data.frame(lon = -149.84, lat = -17.53) #set map coordinates
map1 <- get_map(moorea, zoom = 12, maptype = "satellite") #set map zoom and type

hiSi <- data_na %>%
  filter(Silicate_May >= 1.5, 
         Silicate_August >= 1.5) %>%
  select(c(Site_Number, Lat, Lon, Silicate_May, Silicate_August))

cluster_map <- ggmap(map1) +
  geom_point(data = hiSi,
             aes(x = Lon, y = Lat), color = "white", size = 3) +
  geom_point(data = full_data1, 
             aes(x = Lon, y = Lat, 
                 color = cluster), 
            size = 3) + #set alpha and color aesthetics
  labs(x = "Longitude", y = "Latitude") #label x and y axes

cluster_map

### Create boxplot faceted by clusters
#First, pivot longer
full_data1_longer <- full_data1 %>%
  pivot_longer(cols = c(6:15), 
               names_to = "Nut_parameters", 
               values_to = "Nut_values")

ggplot(full_data1_longer, mapping = aes(x = Nut_parameters,
                                 y= Nut_values + .1, color = cluster)) +
  geom_boxplot()+
  coord_trans(y="log")+
  facet_wrap(~cluster) +
  ylab("Concentration")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = -.1)) 



### Redo everything without variable selection #####################
#Estimated partition
fitted(res_without)

#Estimated probabilities of classification
head(fitted(res_without, type = "probability"))

#Summary
summary(res_without)

#Discriminative power of the variables
plot(res_without)


### Save the clusters as a column in new object

full_data2 <- full_data
full_data2$cluster <- fitted(res_without)

full_data2 <- full_data2 %>%
  mutate("cluster" = factor(cluster))

### Map clusters 
moorea <- data.frame(lon = -149.84, lat = -17.53) #set map coordinates
map1 <- get_map(moorea, zoom = 12, maptype = "satellite") #set map zoom and type

#New dataset with just high silicate values
hiSi <- data_na %>%
  filter(Silicate_May >= 1.5, 
         Silicate_August >= 1.5) %>%
  select(c(Site_Number, Lat, Lon, Silicate_May, Silicate_August))

colors <- c("red", "green", "orange")

cluster_map_without <- ggmap(map1) +
  geom_point(data = hiSi,
             aes(x = Lon, y = Lat), color = "white", size = 3) +
  geom_point(data = full_data2, 
             aes(x = Lon, y = Lat, 
                 color = cluster), 
             size = 3) + #set alpha and color aesthetics
  scale_color_manual(values = colors) +
  labs(x = "Longitude", y = "Latitude", title = "LCM Clustering Map- Jamie") #label x and y axes

cluster_map_without

### Create boxplot faceted by clusters
#First, pivot longer
full_data2_longer <- full_data2 %>%
  pivot_longer(cols = c(6:15), 
               names_to = "Nut_parameters", 
               values_to = "Nut_values")

ggplot(full_data2_longer, mapping = aes(x = Nut_parameters,
                                        y= Nut_values + .1, color = cluster)) +
  geom_boxplot()+
  coord_trans(y="log")+
  facet_wrap(~cluster) +
  ylab("Concentration")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = -.1)) 

