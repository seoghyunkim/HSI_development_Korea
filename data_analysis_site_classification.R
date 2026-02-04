#' # Session information

sessionInfo()
#install.packages("stringr")
#install.packages("knitr")
#install.packages("factoextra")


#' # Load library
rm(list=ls(all.names=T)) # clear all dataddd
pacman::p_load(tidyverse, # load package
               vegan,
               moonBook,
               stringr,
               knitr,
               caret,
               cluster,
               factoextra) 



#' # Read data
#'# Read data 
HSI_DB <- read.csv("./data_outcome/HSI_DB.csv", header=TRUE) %>%
  select(-X) %>% # remove "X" column
  filter(!site == "YEONGSAN_RIVER_MANGWOL") %>% # remove site which is out of streamnetwork
  filter(!site == "YEONGSAN_RIVER_GUJEONG") %>%# remove site which is out of streamnetwork
  filter(!site == "DONGGANG_ST1-2") %>%
  filter(!Channel_ID == "")
#filter(!str_starts(site, "MANCHEON"))

length(unique(HSI_DB$site) ) # site number n = 202


length(unique(HSI_DB$site[HSI_DB$basin == "Eastern" ]))
#'# Data management

#'-Environmental data 
HSI_env <- HSI_DB %>%
  drop_na(substrate2) %>% # remove site which did not measure substrate
  distinct(site_by_season, sample_number, .keep_all = TRUE)%>%# keep single row based on sample number 
  as_tibble()
length(unique(HSI_env$site) ) # site number n = 196 correct





#'- Calculate mean depth 
dep_mean <- aggregate(depth ~ site, data = HSI_env, "mean") %>%
  as_tibble() %>%
  mutate(depth = round(depth, 2)) %>%
  arrange(site)

#'- Calculate mean velocity
vel_mean <- aggregate(velocity~ site, data = HSI_env, "mean") %>%
  as_tibble() %>%
  mutate(velocity = round(velocity, 2)) %>%
  arrange(site)


#'- Calculate percentage of substrate
sub_prop <- HSI_env %>%
  group_by(site) %>%
  count(substrate2) %>% # count 
  mutate(sub_percent = (100 * n/sum(n))) %>% # percentage of each substrate by section
  as_tibble() %>% # change format 
  select(site, substrate2, sub_percent) %>%
  pivot_wider(names_from = substrate2, values_from = sub_percent, values_fill = 0) %>% # section-wdie format and fill 0 value
  mutate(fine_substrate = rowSums(across(c(silt, sand)))) %>%
  mutate(medium_substrate = rowSums(across(c(gravel, pebble)))) %>%
  mutate(coarse_substrate = rowSums(across(c(cobble, boulder)))) %>%
  arrange(site) # change order 

#'- Check result
head(sub_prop, 10)

#' ## Combine depth, velocity, substrate data
env_dat <- dep_mean %>%
  left_join(vel_mean, by ="site") %>%
  left_join(sub_prop, by ="site")

#'- Combine local environment (dep, vel, sub) into site data 
HSI_env2 <- HSI_env %>%
  distinct(site, .keep_all = TRUE) %>%
  left_join(env_dat, by ="site") %>%
  select(basin, site, stream_order, elevation = site_elevation, slope= site_slope, 
         general_slope,  cat_urban= urban, cat_forest= forest,
         cat_agriculture= agriculture , segment_length,ws_area= area_mean,  
         ws_agriculture, ws_forest, ws_urban, 
         depth = depth.y, velocity = velocity.y,
         silt, sand, gravel, pebble, cobble, boulder,
         fine_substrate,medium_substrate = medium_substrate, coarse_substrate)%>%
  mutate(
  cat_agriculture = cat_agriculture * 100,
  cat_forest      = cat_forest * 100,
  cat_urban       = cat_urban * 100,
  ws_forest       = ws_forest * 100,
  ws_urban        = ws_urban * 100,
  ws_agriculture  = ws_agriculture * 100)


length(unique(HSI_env2$site) ) # site number n = 196

summary(HSI_env2)


#' ## Data management for ordination
local_site_env <- HSI_env2 %>%
  select(site, st_order= stream_order, elevation, slope, cat_agriculture, cat_forest, cat_urban,
         ws_area, ws_forest, ws_urban, 
         ws_agriculture, fine_substrate,medium_substrate, coarse_substrate) %>%
  column_to_rownames(var="site")%>%
  drop_na()

attach(local_site_env)  



########## kmeans Clustering#####
#'## (2-1) kmeans clustering

set.seed(123) ##fix the random seed

wss <- sapply(1:10, function(k) {
  kmeans(scale(local_site_env[, 1:13]), centers = k, iter.max = 10000, nstart = 10)$tot.withinss})

elbow_plot <- data.frame(Clusters = 1:10, WSS = wss)
ggplot(elbow_plot, aes(x = Clusters, y = WSS)) +
  geom_line(size = 0.6) +
  geom_point(shape = 16, size = 1.2) +
  labs(title = "Elbow Method", x = "Number of Cluster", y = "WCSS") +
  scale_y_continuous(breaks = seq(0, 2500, 500)) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 10, face = "bold")  
  )


#####Silhouette score#####
packageVersion("cluster")
silhouette_analysis <- function(data, k_values, nstart_values) {
  results <- expand.grid(k = k_values, nstart = nstart_values)  # all combination
  results$silhouette <- NA  # add empty column 
  
  for (i in seq_len(nrow(results))) {
    k <- results$k[i]
    n <- results$nstart[i]
    
    # K-means clustering
    kmeans_result <- kmeans(data, centers = k, nstart = n)
    
    # shihouette calculation
    sil <- silhouette(kmeans_result$cluster, dist(data))
    results$silhouette[i] <- mean(sil[, 3])  # save means of silhoette
  }
  
  return(results)
}

# k value, nstart value setting
k_values <- 2:10          # cluster number
nstart_values <- c(5, 10, 25, 50, 100)  # nstart range

# analysis
results <- silhouette_analysis(local_site_env[, 1:13], k_values, nstart_values)

# result visulization

ggplot(results, aes(x = k, y = silhouette, color = as.factor(nstart), group = nstart)) +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Set1", name = "nstart") +
  labs(
    title = "Silhouette Score",
    x = "Number of Cluster",
    y = "Silhouette Score"
  ) +
  theme_classic()+
  theme(
    axis.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 10, face = "bold")  
  )


optimal_clusters <- 2
site_kmeans <- kmeans(scale(local_site_env[, 1:13]), centers = optimal_clusters, iter.max = 100000)

clustered_data <- data.frame(local_site_env, Cluster = site_kmeans$cluster)


#####Silhouette Width#####
silhouette_vals <- silhouette(site_kmeans$cluster, dist(scale(local_site_env[, 1:12])))
fviz_silhouette(silhouette_vals, palette = c("black", "grey")) +
  labs(title = "Silhouette Plot for K-means Clustering",
       x = "Clusters",
       y = "Silhouette Width") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 3, face = "bold"),  
    axis.text.y = element_text(size = 10, face = "bold"), 
    axis.title = element_text(size = 10, face = "bold")   
  )

silhouette                                                                      

#'# Export final data into data_outcome folder

clustered_grp_data <- clustered_data%>%
  rownames_to_column(var = "site")%>%
  select(site, Cluster)


filename <- paste0("./data_outcome/site_group_dat_k.csv")
write.csv(clustered_grp_data, filename)  



##### PCA_k
env_pca_k <- rda(scale(clustered_data[, 1:13]), scale = TRUE)
summary(env_pca_k)
scores(env_pca_k, display = "sites")   
scores(env_pca_k, display = "species") 

biplot(env_pca_k , type =c("text","point"),xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), scaling= -1)

biplot(env_pca_k , type =c("text","point"),xlim=c(-1.0,1.0), ylim=c(-1.0,1.5), scaling= -1 ,
       col=c("white","black"),
       xlab="PC1 (40.70%)", 
       ylab="PC2 (15.48%)",cex.lab=1.2 )


#'- point
points(env_pca_k, col="blue", select=clustered_data$Cluster=="1", pch=16,cex=0.8, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), scaling= -1)
points(env_pca_k, col="red", select=clustered_data$Cluster=="2", pch=16, cex=0.8,xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), scaling= -1)




ordihull(env_pca_k, clustered_data$Cluster, lty= 1, scaling= -1, lwd= 1, col=c("blue","red"))
legend("bottomleft", legend=c("Group 1","Group 2"),
       col=c("blue","red"), pch=16, bty="n")







#'## Check site group property
group_st_order <- clustered_data %>%
  aggregate(st_order ~Cluster, "range") # variable's range calculation
group_st_order

group_ws_area <- clustered_data %>%
  aggregate(ws_area ~Cluster, "range")
group_ws_area

group_fine_substrate <- clustered_data %>%
  aggregate(fine_substrate ~Cluster, "range")
group_fine_substrate

group_medium_substrate <- clustered_data %>%
  aggregate(medium_substrate ~Cluster, "range")
group_medium_substrate

group_coarse_substrate <- clustered_data %>%
  aggregate(coarse_substrate ~Cluster, "range")
group_coarse_substrate


group_elevation <- clustered_data %>%
  aggregate(elevation ~Cluster, "range")
group_elevation

group_slope <- clustered_data %>%
  aggregate(slope ~Cluster, "range")
group_elevation

group_ws_forest<- clustered_data %>%
  aggregate(ws_forest~Cluster, "range")
group_ws_forest

group_ws_urban<- clustered_data %>%
  aggregate(ws_urban~Cluster, "range")
group_ws_urban

group_ws_agriculture<- clustered_data %>%
  aggregate(ws_agriculture~Cluster, "range")
group_ws_agriculture

group_cat_forest<- clustered_data %>%
  aggregate(cat_forest~Cluster, "range")
group_cat_forest

group_cat_urban<- clustered_data %>%
  aggregate(cat_urban~Cluster, "range")
group_cat_urban


group_cat_agriculture<- clustered_data %>%
  aggregate(cat_agriculture~Cluster, "range")
group_cat_agriculture





#'## Create table

envtable <- mytable(Cluster ~
                      st_order + ws_area+
                      fine_substrate + medium_substrate +coarse_substrate 
                    + slope + elevation 
                    + ws_forest + ws_urban + ws_agriculture 
                    +cat_forest + cat_urban + cat_agriculture 
                    ,data= clustered_data , method= 2, digits= 2)
envtable






